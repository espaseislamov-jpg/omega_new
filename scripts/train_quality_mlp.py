from __future__ import annotations

import argparse
import copy
import json
import math
import random
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import (
    average_precision_score,
    confusion_matrix,
    mean_absolute_error,
    precision_recall_fscore_support,
    precision_recall_curve,
    roc_auc_score,
)
from sklearn.model_selection import GroupKFold

try:
    import torch
    from torch import nn
except ImportError as exc:  # pragma: no cover - exercised in the cloud workflow
    raise SystemExit("PyTorch is required for cloud training: pip install torch") from exc


PROJECT_DIR = Path(__file__).resolve().parents[1]
DEFAULT_DATASET = PROJECT_DIR / "artifacts" / "quality_dataset.csv"
DEFAULT_OUTPUT_DIR = PROJECT_DIR / "artifacts" / "quality_judge"
SEALED_DATES = frozenset({"14072026"})
META_AND_LABEL_COLUMNS = {
    "batch_date",
    "raw_date",
    "sample_id",
    "instrument_no",
    "sample_name",
    "excel_row",
    "match_method",
    "reference",
    "delta",
    "abs_error",
    "error_gt_0_3",
    "error_gt_0_5",
}


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


class FeatureNormalizer:
    def __init__(self) -> None:
        self.medians: np.ndarray | None = None
        self.means: np.ndarray | None = None
        self.scales: np.ndarray | None = None

    def fit(self, values: np.ndarray) -> "FeatureNormalizer":
        values = np.asarray(values, dtype=np.float64)
        with np.errstate(all="ignore"):
            medians = np.nanmedian(values, axis=0)
        medians = np.where(np.isfinite(medians), medians, 0.0)
        filled = np.where(np.isfinite(values), values, medians)
        means = np.mean(filled, axis=0)
        scales = np.std(filled, axis=0)
        scales = np.where(np.isfinite(scales) & (scales > 1e-9), scales, 1.0)
        self.medians = medians.astype(np.float32)
        self.means = means.astype(np.float32)
        self.scales = scales.astype(np.float32)
        return self

    def transform(self, values: np.ndarray) -> np.ndarray:
        if self.medians is None or self.means is None or self.scales is None:
            raise RuntimeError("FeatureNormalizer has not been fitted")
        values = np.asarray(values, dtype=np.float32)
        missing = ~np.isfinite(values)
        filled = np.where(missing, self.medians, values)
        normalized = (filled - self.means) / self.scales
        return np.concatenate([normalized, missing.astype(np.float32)], axis=1).astype(np.float32)


class QualityMLP(nn.Module):
    def __init__(self, input_size: int, hidden_sizes: tuple[int, int] = (64, 32), dropout: float = 0.15) -> None:
        super().__init__()
        self.layer1 = nn.Linear(input_size, hidden_sizes[0])
        self.layer2 = nn.Linear(hidden_sizes[0], hidden_sizes[1])
        self.output = nn.Linear(hidden_sizes[1], 3)
        self.activation = nn.ReLU()
        self.dropout = nn.Dropout(dropout)

    def forward(self, values: torch.Tensor) -> torch.Tensor:
        values = self.dropout(self.activation(self.layer1(values)))
        values = self.dropout(self.activation(self.layer2(values)))
        return self.output(values)


def select_feature_columns(frame: pd.DataFrame) -> list[str]:
    columns: list[str] = []
    for column in frame.columns:
        if column in META_AND_LABEL_COLUMNS:
            continue
        values = pd.to_numeric(frame[column], errors="coerce")
        finite = values[np.isfinite(values)]
        if len(finite) < 2 or finite.nunique(dropna=True) < 2:
            continue
        columns.append(column)
    if not columns:
        raise ValueError("No usable numeric production features found")
    return sorted(columns)


def numeric_matrix(frame: pd.DataFrame, feature_columns: list[str]) -> np.ndarray:
    return np.column_stack([
        pd.to_numeric(frame[column], errors="coerce").to_numpy(dtype=np.float64)
        for column in feature_columns
    ])


def _positive_weight(labels: np.ndarray) -> float:
    positives = float(np.sum(labels > 0.5))
    negatives = float(len(labels) - positives)
    if positives <= 0:
        return 1.0
    return float(np.clip(negatives / positives, 1.0, 25.0))


def multitask_loss(outputs: torch.Tensor, targets: torch.Tensor, weights: tuple[float, float]) -> torch.Tensor:
    regression_loss = nn.functional.smooth_l1_loss(outputs[:, 0], torch.log1p(targets[:, 0]), beta=0.15)
    loss_03 = nn.functional.binary_cross_entropy_with_logits(
        outputs[:, 1], targets[:, 1], pos_weight=outputs.new_tensor(weights[0])
    )
    loss_05 = nn.functional.binary_cross_entropy_with_logits(
        outputs[:, 2], targets[:, 2], pos_weight=outputs.new_tensor(weights[1])
    )
    return regression_loss + 0.55 * loss_03 + 0.90 * loss_05


def train_fold(
    train_x: np.ndarray,
    train_y: np.ndarray,
    valid_x: np.ndarray,
    valid_y: np.ndarray,
    device: torch.device,
    epochs: int,
    patience: int,
    seed: int,
) -> tuple[QualityMLP, int]:
    set_seed(seed)
    model = QualityMLP(train_x.shape[1]).to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=1.5e-3, weight_decay=2e-4)
    train_values = torch.as_tensor(train_x, dtype=torch.float32, device=device)
    train_targets = torch.as_tensor(train_y, dtype=torch.float32, device=device)
    valid_values = torch.as_tensor(valid_x, dtype=torch.float32, device=device)
    valid_targets = torch.as_tensor(valid_y, dtype=torch.float32, device=device)
    class_weights = (_positive_weight(train_y[:, 1]), _positive_weight(train_y[:, 2]))

    best_loss = math.inf
    best_epoch = 1
    best_state = copy.deepcopy(model.state_dict())
    stale_epochs = 0
    for epoch in range(1, epochs + 1):
        model.train()
        optimizer.zero_grad(set_to_none=True)
        loss = multitask_loss(model(train_values), train_targets, class_weights)
        loss.backward()
        nn.utils.clip_grad_norm_(model.parameters(), 5.0)
        optimizer.step()

        model.eval()
        with torch.no_grad():
            valid_loss = float(multitask_loss(model(valid_values), valid_targets, class_weights).cpu())
        if valid_loss < best_loss - 1e-5:
            best_loss = valid_loss
            best_epoch = epoch
            best_state = copy.deepcopy(model.state_dict())
            stale_epochs = 0
        else:
            stale_epochs += 1
            if stale_epochs >= patience:
                break

    model.load_state_dict(best_state)
    return model, best_epoch


def predict(model: QualityMLP, values: np.ndarray, device: torch.device) -> np.ndarray:
    model.eval()
    with torch.no_grad():
        outputs = model(torch.as_tensor(values, dtype=torch.float32, device=device)).cpu().numpy()
    prediction = np.empty_like(outputs, dtype=np.float64)
    prediction[:, 0] = np.expm1(outputs[:, 0]).clip(min=0.0)
    prediction[:, 1:] = 1.0 / (1.0 + np.exp(-outputs[:, 1:]))
    return prediction


def choose_f2_threshold(labels: np.ndarray, probabilities: np.ndarray) -> float:
    if len(np.unique(labels)) < 2:
        return 0.5
    precision, recall, thresholds = precision_recall_curve(labels, probabilities)
    if len(thresholds) == 0:
        return 0.5
    precision = precision[:-1]
    recall = recall[:-1]
    scores = 5.0 * precision * recall / np.maximum(4.0 * precision + recall, 1e-12)
    return float(thresholds[int(np.nanargmax(scores))])


def classification_metrics(labels: np.ndarray, probabilities: np.ndarray, threshold: float) -> dict:
    predicted = (probabilities >= threshold).astype(int)
    precision, recall, f1, _ = precision_recall_fscore_support(
        labels, predicted, average="binary", zero_division=0
    )
    matrix = confusion_matrix(labels, predicted, labels=[0, 1])
    result = {
        "threshold": float(threshold),
        "precision": float(precision),
        "recall": float(recall),
        "f1": float(f1),
        "confusion_matrix": matrix.tolist(),
        "positives": int(np.sum(labels)),
    }
    if len(np.unique(labels)) >= 2:
        result["roc_auc"] = float(roc_auc_score(labels, probabilities))
        result["average_precision"] = float(average_precision_score(labels, probabilities))
    return result


def export_numpy_model(
    path: Path,
    model: QualityMLP,
    normalizer: FeatureNormalizer,
    feature_columns: list[str],
    threshold_03: float,
    threshold_05: float,
) -> None:
    state = {key: value.detach().cpu().numpy() for key, value in model.state_dict().items()}
    np.savez_compressed(
        path,
        feature_names=np.asarray(feature_columns, dtype=str),
        medians=normalizer.medians,
        means=normalizer.means,
        scales=normalizer.scales,
        threshold_03=np.asarray([threshold_03], dtype=np.float32),
        threshold_05=np.asarray([threshold_05], dtype=np.float32),
        **{f"state_{key.replace('.', '__')}": value for key, value in state.items()},
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Train the cloud-only multitask neural quality judge.")
    parser.add_argument("--dataset", type=Path, default=DEFAULT_DATASET)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--epochs", type=int, default=600)
    parser.add_argument("--patience", type=int, default=70)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--device", choices=["auto", "cpu", "cuda"], default="auto")
    args = parser.parse_args()

    # Keep accidental local execution gentle; cloud GPU use is unaffected.
    torch.set_num_threads(min(2, max(1, torch.get_num_threads())))
    frame = pd.read_csv(args.dataset)
    if frame.empty:
        raise SystemExit("Training dataset is empty")
    dates = set(frame["batch_date"].astype(str))
    leaked = dates.intersection(SEALED_DATES)
    if leaked:
        raise SystemExit(f"Sealed batches found in training data: {sorted(leaked)}")
    if frame["batch_date"].nunique() < 3:
        raise SystemExit("At least three complete batches are required for grouped validation")

    device_name = "cuda" if args.device == "auto" and torch.cuda.is_available() else args.device
    if device_name == "auto":
        device_name = "cpu"
    device = torch.device(device_name)
    feature_columns = select_feature_columns(frame)
    raw_x = numeric_matrix(frame, feature_columns)
    y = np.column_stack([
        pd.to_numeric(frame["abs_error"], errors="raise").to_numpy(dtype=np.float32),
        pd.to_numeric(frame["error_gt_0_3"], errors="raise").to_numpy(dtype=np.float32),
        pd.to_numeric(frame["error_gt_0_5"], errors="raise").to_numpy(dtype=np.float32),
    ])
    groups = frame["batch_date"].astype(str).to_numpy()

    splits = min(5, len(np.unique(groups)))
    splitter = GroupKFold(n_splits=splits)
    oof = np.full((len(frame), 3), np.nan, dtype=np.float64)
    best_epochs: list[int] = []
    fold_rows: list[dict] = []
    for fold, (train_idx, valid_idx) in enumerate(splitter.split(raw_x, y[:, 2], groups), start=1):
        normalizer = FeatureNormalizer().fit(raw_x[train_idx])
        train_x = normalizer.transform(raw_x[train_idx])
        valid_x = normalizer.transform(raw_x[valid_idx])
        model, best_epoch = train_fold(
            train_x,
            y[train_idx],
            valid_x,
            y[valid_idx],
            device=device,
            epochs=args.epochs,
            patience=args.patience,
            seed=args.seed + fold,
        )
        oof[valid_idx] = predict(model, valid_x, device)
        best_epochs.append(best_epoch)
        fold_rows.append({
            "fold": fold,
            "train_rows": int(len(train_idx)),
            "valid_rows": int(len(valid_idx)),
            "valid_batches": sorted(set(groups[valid_idx])),
            "best_epoch": int(best_epoch),
        })
        print(f"fold={fold}/{splits} best_epoch={best_epoch} valid_batches={fold_rows[-1]['valid_batches']}", flush=True)

    if not np.isfinite(oof).all():
        raise RuntimeError("OOF predictions are incomplete")
    threshold_03 = choose_f2_threshold(y[:, 1].astype(int), oof[:, 1])
    threshold_05 = choose_f2_threshold(y[:, 2].astype(int), oof[:, 2])
    report = {
        "rows": int(len(frame)),
        "batches": sorted(dates),
        "sealed_dates_excluded": sorted(SEALED_DATES),
        "feature_count": len(feature_columns),
        "device": str(device),
        "grouped_folds": fold_rows,
        "oof_expected_error_mae": float(mean_absolute_error(y[:, 0], oof[:, 0])),
        "oof_error_gt_0_3": classification_metrics(y[:, 1].astype(int), oof[:, 1], threshold_03),
        "oof_error_gt_0_5": classification_metrics(y[:, 2].astype(int), oof[:, 2], threshold_05),
    }

    final_epochs = max(20, int(np.median(best_epochs)))
    final_normalizer = FeatureNormalizer().fit(raw_x)
    final_x = final_normalizer.transform(raw_x)
    set_seed(args.seed)
    final_model = QualityMLP(final_x.shape[1]).to(device)
    optimizer = torch.optim.AdamW(final_model.parameters(), lr=1.5e-3, weight_decay=2e-4)
    values = torch.as_tensor(final_x, dtype=torch.float32, device=device)
    targets = torch.as_tensor(y, dtype=torch.float32, device=device)
    class_weights = (_positive_weight(y[:, 1]), _positive_weight(y[:, 2]))
    for _ in range(final_epochs):
        final_model.train()
        optimizer.zero_grad(set_to_none=True)
        loss = multitask_loss(final_model(values), targets, class_weights)
        loss.backward()
        nn.utils.clip_grad_norm_(final_model.parameters(), 5.0)
        optimizer.step()
    report["final_epochs"] = final_epochs

    args.output_dir.mkdir(parents=True, exist_ok=True)
    export_numpy_model(
        args.output_dir / "quality_judge_model.npz",
        final_model,
        final_normalizer,
        feature_columns,
        threshold_03,
        threshold_05,
    )
    frame.assign(
        predicted_abs_error_oof=oof[:, 0],
        probability_gt_0_3_oof=oof[:, 1],
        probability_gt_0_5_oof=oof[:, 2],
    ).to_csv(args.output_dir / "quality_judge_oof_predictions.csv", index=False)
    (args.output_dir / "quality_judge_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False),
        encoding="utf-8",
    )
    print(json.dumps(report, ensure_ascii=False, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
