"""Tkinter GUI for the production Omega chromatogram engine."""

import json
import re
import logging
import sys
from concurrent.futures import Future
from threading import Thread
from pathlib import Path

from omega_path_compat import configure_windows_path_compat

configure_windows_path_compat()

import numpy as np
import pandas as pd
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk

import omega_core
from omega_core import metrics as core_metrics
from omega_core import instrument_profiles
from omega_core import rt_profile, manual_edits, runtime
from omega_core.results import finalize_result
from omega_core.io import ensure_runtime_file
from omega_ui import WorkspaceUI, readable_status
from omega_version import APP_NAME
from omega_workflow import WorkflowUI
from omega_profiles_ui import ProfileManagerUI
from omega_manual_ui import ManualIntegrationUI
from omega_plot_ui import ChromatogramPlotUI
# Retain the historical helper import for external diagnostic tools.
from omega_plot_ui import _get_x_column_name as _get_x_column_name


def process_chromatogram_batch(dataframe: pd.DataFrame, reference_targets: pd.DataFrame) -> dict:
    """Run the same modular engine used by the regression harness."""
    engine_input = dataframe.copy()
    if "x_corrected" not in engine_input.columns and "x" in engine_input.columns:
        engine_input["x_corrected"] = pd.to_numeric(engine_input["x"], errors="coerce")
    result = dict(omega_core.process_batch(engine_input, reference_targets))
    result["engine"] = "omega_core"
    result.setdefault("total_area", result.get("omega", {}).get("total_area", np.nan))
    return result


class ChromatogramToolbar(NavigationToolbar2Tk):
    """Keep viewing/export tools; subplot layout is an implementation detail."""
    toolitems = tuple(item for item in NavigationToolbar2Tk.toolitems
                      if item[3] != "configure_subplots")


class ChromatogramApp(WorkflowUI, WorkspaceUI, ProfileManagerUI, ManualIntegrationUI, ChromatogramPlotUI):
    def __init__(self, root: tk.Tk):
        self.root = root
        self.root.report_callback_exception = self.report_callback_exception
        self.root.title(f"{APP_NAME} — анализ хроматограмм")
        screen_width = max(int(self.root.winfo_screenwidth()), 800)
        screen_height = max(int(self.root.winfo_screenheight()), 600)
        self.screen_width = screen_width
        self.screen_height = screen_height
        self.compact_ui = screen_width < 1650 or screen_height < 950
        initial_width = min(1780, max(960, int(screen_width * 0.94)))
        initial_height = min(1040, max(640, int(screen_height * 0.88)))
        initial_width = min(initial_width, screen_width)
        initial_height = min(initial_height, screen_height)
        self.initial_window_width = initial_width
        self.initial_window_height = initial_height
        self.root.geometry(f"{initial_width}x{initial_height}")
        min_width = max(900, min(1050, int(screen_width * 0.72)))
        min_height = max(580, min(680, int(screen_height * 0.70)))
        self.root.minsize(min_width, min_height)

        self.reference_json_path = ensure_runtime_file("reference_targets_reverted_c22fixed.json")
        self.base_reference_targets = omega_core.load_reference_targets(self.reference_json_path)
        self.profile_store = instrument_profiles.load_store(self.base_reference_targets)
        self.active_instrument_profile = instrument_profiles.active_profile(self.profile_store)
        self.reference_targets = instrument_profiles.apply_profile_to_targets(
            self.base_reference_targets, self.active_instrument_profile
        )

        self.current_file = None
        self.current_sample_name = ""
        self.loaded_batches = []
        self.current_batch_index = 0
        self.batch_tree = None
        self.batch_results_window = None
        self.batch_results_tree = None
        self._batch_tree_syncing = False
        self._preload_batch_index = 0
        self._preload_after_id = None
        self._background_future = None
        self.batch_progress_window = None
        self.batch_progress_label_var = tk.StringVar(value="")
        self.batch_progress_detail_var = tk.StringVar(value="")
        self.batch_progress_bar = None
        self.df_processed = None
        self.best_window = None
        self.peaks_df = pd.DataFrame()
        self.matched_targets_df = pd.DataFrame()
        self.current_rt_shift = 0.0
        self.selected_target_code = None
        self.manual_start_var = tk.StringVar(value="")
        self.manual_end_var = tk.StringVar(value="")
        self._manual_drag_active_boundary = None
        self._manual_drag_pending_bounds = None
        self._manual_drag_after_id = None
        self._manual_drag_axis = None
        self._manual_pick = None
        self._manual_pick_artist = None
        self._manual_overlay_artists = {}

        self.status_var = tk.StringVar(value="Выбери CSV-файл.")
        self.file_var = tk.StringVar(value="Файл не выбран")
        self.omega_var = tk.StringVar(value="—")
        self.integration_var = tk.StringVar(value="Integration: —")
        self.gamma_var = tk.StringVar(value="γ-Linolenic: —")
        self.batch_var = tk.StringVar(value="Проба: —")
        self.confidence_var = tk.StringVar(value="Качество пиков: —")
        self.current_confidence = None
        self.profile_var = tk.StringVar(value=self.active_instrument_profile["name"])
        self.profile_combo = None
        self.profile_window = None

        self._build_ui()
        self.setup_workflow()

    def _build_ui(self):
        self.build_workspace(ChromatogramToolbar)

    def export_boundary_history(self):
        if not self.loaded_batches:
            messagebox.showinfo("История границ", "Сначала откройте CSV.", parent=self.root)
            return
        batch = self.loaded_batches[self.current_batch_index]
        path = filedialog.asksaveasfilename(
            parent=self.root, title="Сохранить историю автоматических границ",
            defaultextension=".json", initialfile="omega_boundary_history.json",
            filetypes=[("JSON", "*.json")],
        )
        if not path:
            return
        payload = {
            "scope": "automatic_cluster_refinement_before_manual_edits",
            "sample": self.current_sample_name,
            "profile": self.active_instrument_profile,
            "baseline_mode": batch.get("baseline_mode"),
            "changes": batch.get("boundary_history", []),
            "manual_changes": batch.get("manual_history", []),
        }
        try:
            with open(path, "w", encoding="utf-8") as stream:
                json.dump(payload, stream, ensure_ascii=False, indent=2, allow_nan=False)
        except (OSError, TypeError, ValueError) as exc:
            messagebox.showerror("История границ", str(exc), parent=self.root)
            return
        self.status_var.set("История автоматических границ сохранена.")

    def _apply_responsive_pane_positions(self):
        try:
            self.position_workspace()
        except tk.TclError:
            pass

    def update_batch_navigation(self):
        total = len(self.loaded_batches)
        if total <= 0:
            self.batch_var.set("Проба: —")
            self.prev_button.state(["disabled"])
            self.next_button.state(["disabled"])
            self.batch_results_button.state(["disabled"])
            self.confidence_var.set("Уверенность: —")
            self.current_confidence = None
            if self.confidence_button is not None:
                self.confidence_button.state(["disabled"])
            if self.batch_tree is not None:
                for item_id in self.batch_tree.get_children():
                    self.batch_tree.delete(item_id)
            return

        self.batch_var.set(f"Проба: {self.current_batch_index + 1}/{total}")
        self.batch_results_button.state(["!disabled"])
        if self.current_batch_index > 0:
            self.prev_button.state(["!disabled"])
        else:
            self.prev_button.state(["disabled"])

        if self.current_batch_index < total - 1:
            self.next_button.state(["!disabled"])
        else:
            self.next_button.state(["disabled"])
        self.populate_main_batch_tree()

    def show_confidence_details(self):
        confidence = self.current_confidence
        if not confidence or not np.isfinite(confidence.get("score", np.nan)):
            messagebox.showinfo("Качество пиков", "Нет данных для оценки геометрии пиков.", parent=self.root)
            return

        risk = confidence.get("high_error_risk", {})
        risk_score = risk.get("score", 0) if isinstance(risk, dict) else 0
        geometry_score = pd.to_numeric(
            pd.Series([confidence.get("geometry_score")]), errors="coerce"
        ).iloc[0]
        emergency_judge = bool(
            risk.get("instrument_profile_judge_emergency_enabled", False)
        ) if isinstance(risk, dict) else False
        if risk_score >= 95:
            recommendation = "СТОП: откройте отмеченные пики и переинтегрируйте их вручную"
        elif risk_score >= 85:
            recommendation = (
                "ПРОВЕРИТЬ: экстренный судья нашёл риск крупной ошибки; "
                "сверьте отмеченные пики"
                if emergency_judge
                else "ПРОВЕРИТЬ: откройте отмеченные пики; при неверной границе переинтегрируйте"
            )
        elif np.isfinite(geometry_score) and geometry_score < core_metrics.GEOMETRY_STOP_SCORE:
            recommendation = "СТОП: геометрия пиков требует ручной проверки"
        elif np.isfinite(geometry_score) and geometry_score < core_metrics.GEOMETRY_READY_SCORE:
            recommendation = "ПРОВЕРИТЬ: геометрия пиков недостаточно уверенная"
        else:
            recommendation = (
                "ГОТОВО: экстренная проверка не нашла грубых признаков ошибки"
                if emergency_judge
                else "Проверки геометрии пройдены; это не подтверждение аналитической точности"
            )
        lines = [
            f"Рекомендация: {recommendation}",
            "",
        ]
        reasons = confidence.get("reasons") or []
        if reasons:
            lines.append("На что обратить внимание:")
            lines.extend(reasons)
        else:
            lines.append("Заметных проблем не найдено.")

        risky_codes = risk.get("peak_codes", []) if isinstance(risk, dict) else []
        if risky_codes and not self.matched_targets_df.empty:
            lines.extend(["", "Площади пиков для проверки:"])
            for code in risky_codes:
                target = self.matched_targets_df[self.matched_targets_df["code"] == code]
                if target.empty:
                    continue
                area = pd.to_numeric(target["area"], errors="coerce").iloc[0]
                area_text = f"{float(area):,.1f}".replace(",", " ") if np.isfinite(area) else "—"
                lines.append(f"• {code}: {area_text}")

        messagebox.showinfo("Качество пиков", "\n".join(lines), parent=self.root)


    def export_runtime_diagnostics(self):
        path = filedialog.asksaveasfilename(parent=self.root, defaultextension=".json", initialfile="omega_diagnostics.json")
        if not path:
            return
        payload = runtime.snapshot(self.reference_targets)
        payload["profile"] = self.active_instrument_profile
        payload["input_sha256"] = manual_edits.file_digest(self.current_file) if self.current_file else None
        payload["samples"] = [{"sample": b["sample_name"], "baseline": b.get("baseline_mode"),
                               "report_status": b.get("report_status"), "reasons": b.get("report_reasons", [])}
                              for b in self.loaded_batches]
        Path(path).write_text(json.dumps(payload, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")
        self.status_var.set("Диагностика окружения и расчёта сохранена.")


    def handle_batch_tree_selection(self, event=None):
        if self._batch_tree_syncing or self.batch_tree is None:
            return
        selection = self.batch_tree.selection()
        if not selection:
            return
        target_index = int(selection[0])
        if target_index != self.current_batch_index:
            self.load_batch_at_index(target_index)

    def prev_batch(self):
        visible = [int(i) for i in self.batch_tree.get_children() if int(i) < self.current_batch_index]
        if visible:
            self.load_batch_at_index(visible[-1])

    def next_batch(self):
        visible = [int(i) for i in self.batch_tree.get_children() if int(i) > self.current_batch_index]
        if visible:
            self.load_batch_at_index(visible[0])

    def process_batch(self, batch: dict):
        if batch.get("processed_df") is not None:
            return batch
        batch.update(process_chromatogram_batch(batch["dataframe"], self.reference_targets))
        return batch

    def load_batch_at_index(self, index: int):
        self.cancel_manual_pick()
        if index < 0 or index >= len(self.loaded_batches):
            return

        preserve_plot_view = self.df_processed is not None
        batch = self.loaded_batches[index]
        self.process_batch(batch)
        self.current_batch_index = index
        self.current_sample_name = batch["sample_name"]
        self.file_var.set(
            f"{self.current_file.name}  /  {self.current_sample_name}"
        )

        self.sample_title_var.set(f"Проба {self.current_sample_name}")
        self.df_processed = batch["processed_df"]
        self.best_window = batch["best_window"]
        self.peaks_df = batch["peaks_df"]
        self.matched_targets_df = batch["matched_targets_df"]
        self.current_rt_shift = batch["rt_shift"]
        self.preview_specs = rt_profile.preview_windows(self.matched_targets_df)
        available_codes = set(self.matched_targets_df.get("code", pd.Series(dtype=str)).astype(str))
        if self.selected_target_code in available_codes:
            self.load_selected_integration_bounds(silent=True)
        else:
            self.selected_target_code = None
            self.manual_start_var.set("")
            self.manual_end_var.set("")
        self.update_batch_navigation()
        self.refresh_peaks(preserve_plot_view=preserve_plot_view)

    @staticmethod
    def _configure_quality_tags(tree: ttk.Treeview):
        tree.tag_configure("quality_stop", background="#FCE4E7", foreground="#A52234")
        tree.tag_configure("quality_check", background="#FFF5DF", foreground="#805500")
        tree.tag_configure("quality_good", background="#FFFFFF", foreground="#25714A")
        tree.tag_configure("quality_pending", foreground="#6f7780")

    @staticmethod
    def _operator_quality(confidence: dict | None) -> tuple[str, str]:
        """Translate the judge into one unambiguous operator action."""
        if not isinstance(confidence, dict) or not confidence:
            return "ОЖИДАНИЕ — ещё не рассчитано", "quality_pending"
        risk = confidence.get("high_error_risk", {})
        risk_score = risk.get("score", 0) if isinstance(risk, dict) else 0
        geometry_score = pd.to_numeric(
            pd.Series([confidence.get("geometry_score")]), errors="coerce"
        ).iloc[0]
        retry = confidence.get("structural_retry", {})
        if risk_score >= 95:
            return "СТОП — переинтегрировать", "quality_stop"
        if risk_score >= 85:
            reason_codes = set(risk.get("reason_codes", [])) if isinstance(risk, dict) else set()
            if "instrument_profile_uncalibrated" in reason_codes:
                return "ПРОВЕРИТЬ — судья профиля ещё не настроен", "quality_check"
            if "emergency_profile_geometry" in reason_codes:
                risky_codes = ", ".join(risk.get("peak_codes", [])) if isinstance(risk, dict) else ""
                return f"ПРОВЕРИТЬ — экстренно: {risky_codes or 'границы'}", "quality_check"
            if risk.get("review_priority") == "c22_first":
                return "ПРОВЕРИТЬ СНАЧАЛА — C22:5, C22:4", "quality_check"
            risky_codes = ", ".join(risk.get("peak_codes", [])) if isinstance(risk, dict) else ""
            suffix = f": {risky_codes}" if risky_codes else " отмеченные пики"
            return f"ПРОВЕРИТЬ{suffix}", "quality_check"
        if np.isfinite(geometry_score) and geometry_score < core_metrics.GEOMETRY_STOP_SCORE:
            return "СТОП — ручная проверка пиков", "quality_stop"
        if np.isfinite(geometry_score) and geometry_score < core_metrics.GEOMETRY_READY_SCORE:
            return "ПРОВЕРИТЬ — геометрию пиков", "quality_check"
        if isinstance(retry, dict) and retry.get("accepted"):
            return "ГОТОВО — перепроверено автоматически", "quality_good"
        if isinstance(risk, dict) and risk.get("instrument_profile_judge_emergency_enabled"):
            return "ГОТОВО — экстренная проверка пройдена", "quality_good"
        return "ПРОВЕРКИ ПРОЙДЕНЫ", "quality_good"

    def build_batch_results_rows(self, process_all: bool = False):
        def percent_text(batch: dict, code: str) -> str:
            matched = batch.get("matched_targets_df")
            if not isinstance(matched, pd.DataFrame) or matched.empty:
                return ""
            row = matched.loc[matched["code"].astype(str) == code]
            if row.empty:
                return ""
            value = pd.to_numeric(row["percent_area"], errors="coerce").iloc[0]
            return f"{float(value):.4f}" if np.isfinite(value) else ""

        rows = []
        for index, batch in enumerate(self.loaded_batches):
            if process_all:
                self.process_batch(batch)
            omega_report = batch.get("omega_report")
            if omega_report is None and isinstance(batch.get("omega"), dict):
                omega_report = batch.get("omega", {}).get("omega3_trio", np.nan)
            value = omega_report if omega_report is not None else np.nan
            value_text = f"{value:.4f}" if np.isfinite(value) else ""
            confidence = batch.get("confidence") if isinstance(batch.get("confidence"), dict) else {}
            confidence_text, quality_tag = self._operator_quality(confidence)
            rows.append((
                index,
                batch.get("sample_name", f"Batch {index + 1}"),
                value_text,
                confidence_text,
                percent_text(batch, "C20:4N6"),
                percent_text(batch, "C18:2N6C"),
                percent_text(batch, "C18:3N6"),
                quality_tag,
            ))
        return rows

    @staticmethod
    def _sample_number(sample_name: str) -> str:
        text = str(sample_name or "").strip()
        if re.match(r"^O\d+_", text, flags=re.IGNORECASE):
            text = text.split("_", 1)[1]
        if text.upper().endswith(".D"):
            text = text[:-2]
        return text

    def _populate_batch_tree_widget(self, tree: ttk.Treeview, process_all: bool = False):
        if tree is None:
            return
        selected_iid = str(self.current_batch_index) if self.loaded_batches else None
        for item_id in tree.get_children():
            tree.delete(item_id)
        show_confidence = "confidence" in set(tree["columns"])
        for (
            index,
            sample_name,
            value_text,
            confidence_text,
            arachidonic_text,
            linoleic_text,
            gamma_linolenic_text,
            quality_tag,
        ) in self.build_batch_results_rows(process_all=process_all):
            if tree is self.batch_tree:
                query = self.sample_search_var.get().strip().casefold()
                if query and query not in str(sample_name).casefold():
                    continue
                if self.review_only_var.get() and quality_tag not in ('quality_check','quality_stop'):
                    continue
                selected_filter = self.review_filter_var.get()
                if selected_filter == 'Требуют проверки' and quality_tag not in ('quality_check','quality_stop'):
                    continue
                if selected_filter == 'Требуют обязательной проверки' and quality_tag != 'quality_stop':
                    continue
                confidence_text = {'quality_good':'Готово', 'quality_check':'Проверить', 'quality_stop':'Стоп', 'quality_pending':'Ожидание'}[quality_tag]
                display_name = str(sample_name).split('_', 1)[0] if re.match(r'^O\d+_', str(sample_name), flags=re.IGNORECASE) else str(sample_name)
                value_text = f'{float(value_text):.2f}' if value_text else '?'
            else:
                display_name = self._sample_number(sample_name)
            values = (
                display_name,
                value_text,
                confidence_text,
                arachidonic_text,
                linoleic_text,
                gamma_linolenic_text,
            ) if show_confidence else (display_name, value_text)
            if tree is self.batch_results_tree and show_confidence:
                values = (display_name, value_text, arachidonic_text,
                          linoleic_text, gamma_linolenic_text, confidence_text)
            tree.insert("", "end", iid=str(index), values=values, tags=(quality_tag,))
        if selected_iid is not None and tree.exists(selected_iid):
            self._batch_tree_syncing = True
            tree.selection_set(selected_iid)
            tree.focus(selected_iid)
            tree.see(selected_iid)
            self._batch_tree_syncing = False

    def populate_main_batch_tree(self):
        self._populate_batch_tree_widget(self.batch_tree, process_all=False)
        if hasattr(self, 'sample_count_var'):
            count = len(self.batch_tree.get_children())
            self.sample_count_var.set(f'Показано {count} из {len(self.loaded_batches)}')
            self.filter_description.configure(text='Красные пробы: обязательная проверка' if self.review_filter_var.get() == 'Требуют обязательной проверки' else '')
            visible = [int(i) for i in self.batch_tree.get_children()]
            self.prev_button.state(['!disabled'] if any(i < self.current_batch_index for i in visible) else ['disabled'])
            self.next_button.state(['!disabled'] if any(i > self.current_batch_index for i in visible) else ['disabled'])

    def copy_batch_results(self, selected_only: bool):
        if self.batch_results_tree is None:
            return

        item_ids = list(self.batch_results_tree.selection()) if selected_only else list(self.batch_results_tree.get_children())
        if not item_ids:
            return

        lines = []
        for item_id in item_ids:
            values = self.batch_results_tree.item(item_id, "values")
            lines.append("\t".join(str(v) for v in values))
        self.root.clipboard_clear()
        self.root.clipboard_append("\n".join(lines))
        self.status_var.set(f"Скопировано строк: {len(lines)}")

    def handle_batch_results_copy(self, event):
        if self.batch_results_tree is None:
            return
        self.copy_batch_results(selected_only=True)

    def jump_to_batch_from_results(self, event=None):
        if self.batch_results_tree is None:
            return
        selection = self.batch_results_tree.selection()
        if not selection:
            return
        target_index = int(selection[0])
        self.load_batch_at_index(target_index)

    def populate_batch_results_tree(self):
        # Do not process the whole file from a Tk callback.  A large field batch
        # can take many seconds and blocks every window event while it runs.
        self._populate_batch_tree_widget(self.batch_results_tree, process_all=False)

    def open_batch_results_window(self):
        if not self.loaded_batches:
            return

        if self.batch_results_window is not None and self.batch_results_window.winfo_exists():
            self.populate_batch_results_tree()
            self.batch_results_window.deiconify()
            self.batch_results_window.lift()
            self.batch_results_window.focus_force()
            return

        self.batch_results_window = tk.Toplevel(self.root)
        self.batch_results_window.title("Результаты серии")
        self.batch_results_window.geometry("1200x650")
        self.batch_results_window.minsize(900, 420)

        frame = ttk.Frame(self.batch_results_window, padding=10)
        frame.pack(fill="both", expand=True)

        toolbar = ttk.Frame(frame)
        toolbar.pack(fill="x", pady=(0, 8))
        ttk.Button(toolbar, text="Копировать выбранные", command=lambda: self.copy_batch_results(selected_only=True)).pack(side="left")
        ttk.Button(toolbar, text="Копировать всё", command=lambda: self.copy_batch_results(selected_only=False)).pack(side="left", padx=(8, 0))
        ttk.Button(toolbar, text="Открыть пробу", command=self.jump_to_batch_from_results).pack(side="left", padx=(8, 0))

        tree_frame = ttk.Frame(frame)
        tree_frame.pack(fill="both", expand=True)
        columns = (
            "sample_name",
            "omega_value",
            "arachidonic",
            "linoleic",
            "gamma_linolenic",
            "confidence",
        )
        self.batch_results_tree = ttk.Treeview(tree_frame, columns=columns, show="headings", height=18, selectmode="extended")
        self.batch_results_tree.heading("sample_name", text="Номер образца")
        self.batch_results_tree.heading("omega_value", text="Значение")
        self.batch_results_tree.heading("confidence", text="Проверка")
        self.batch_results_tree.heading("arachidonic", text="Арахидоновая, %")
        self.batch_results_tree.heading("linoleic", text="Линолевая, %")
        self.batch_results_tree.heading("gamma_linolenic", text="γ-Линоленовая, %")
        self.batch_results_tree.column("sample_name", width=220, minwidth=130, anchor="w", stretch=True)
        self.batch_results_tree.column("omega_value", width=110, minwidth=80, anchor="center", stretch=False)
        self.batch_results_tree.column("confidence", width=300, minwidth=220, anchor="w", stretch=True)
        self.batch_results_tree.column("arachidonic", width=125, minwidth=105, anchor="center", stretch=False)
        self.batch_results_tree.column("linoleic", width=115, minwidth=95, anchor="center", stretch=False)
        self.batch_results_tree.column("gamma_linolenic", width=135, minwidth=115, anchor="center", stretch=False)
        self.batch_results_tree.pack(side="left", fill="both", expand=True)
        self._configure_quality_tags(self.batch_results_tree)

        scrollbar = ttk.Scrollbar(tree_frame, orient="vertical", command=self.batch_results_tree.yview)
        scrollbar.pack(side="right", fill="y")
        self.batch_results_tree.configure(yscrollcommand=scrollbar.set)
        self.batch_results_tree.bind("<Double-1>", self.jump_to_batch_from_results)
        self.batch_results_tree.bind("<Control-c>", self.handle_batch_results_copy)

        help_label = ttk.Label(frame, text="Ctrl+C копирует выделенные строки. Двойной клик открывает выбранный образец.")
        help_label.pack(fill="x", pady=(8, 0))

        def on_close():
            if self.batch_results_window is not None:
                self.batch_results_window.destroy()
            self.batch_results_window = None
            self.batch_results_tree = None

        self.batch_results_window.protocol("WM_DELETE_WINDOW", on_close)
        self.populate_batch_results_tree()


    def open_file(self, profile_id=None):
        if self._background_future is not None or self._preload_after_id is not None:
            return
        try:
            file_path = filedialog.askopenfilename(parent=self.root, title="Выберите CSV", filetypes=[("CSV", "*.csv *.CSV"), ("All", "*.*")])
            if not file_path:
                return
            logging.getLogger("omega").info("Reading CSV: %s", file_path)
            self.show_batch_progress_window()
            self.batch_progress_label_var.set("Чтение CSV…")
            self.batch_progress_detail_var.set(Path(file_path).name)
            self.batch_progress_bar.configure(mode="indeterminate")
            self.batch_progress_bar.start(50)
            self.status_var.set("Чтение CSV…")
            self._run_background(
                lambda: omega_core.load_batches(Path(file_path), cutoff_minutes=4.0),
                lambda batches: self.finish_new_calculation(Path(file_path), batches, profile_id),
            )
        except Exception:
            self.report_callback_exception(*sys.exc_info())

    def _finish_file_load(self, file_path, batches):
        if not batches:
            raise ValueError("В CSV нет проб для анализа.")
        self.current_file = file_path
        self.loaded_batches = batches
        # A newly read file must never retain a previous file's plot or number.
        self.df_processed = None
        self.peaks_df = pd.DataFrame()
        self.matched_targets_df = pd.DataFrame()
        self.current_confidence = None
        self.selected_target_code = None
        self.current_sample_name = ""
        self.current_batch_index = 0
        self.omega_var.set("Расчёт…")
        self.manual_start_var.set("")
        self.manual_end_var.set("")
        self.sample_title_var.set("Загрузка проб…")
        self.sample_search_var.set("")
        self.review_only_var.set(False)
        self.review_filter_var.set('Все пробы')
        self.update_review_panel()
        self.confidence_button.state(['disabled'])
        self.file_var.set(f"Файл: {file_path.name}")
        for tree in (self.tree, self.batch_tree, self.batch_results_tree):
            if tree is not None and tree.winfo_exists():
                tree.delete(*tree.get_children())
        for axis in self.figure.axes:
            axis.clear()
        self._manual_overlay_artists = {}
        self.canvas.draw_idle()
        logging.getLogger("omega").info("CSV read complete: %d samples", len(batches))
        self._preload_batch_index = 0
        self.batch_progress_bar.stop()
        self.batch_progress_bar.configure(mode="determinate", maximum=len(batches), value=0)
        self.status_var.set(f"Загружено проб: {len(batches)}. Запускаю расчёт.")
        self._preload_after_id = self.root.after(25, self.preload_loaded_batches)

    def report_callback_exception(self, exc_type, exc, tb):
        logging.getLogger("omega").error("GUI operation failed", exc_info=(exc_type, exc, tb))
        self.close_batch_progress_window()
        self.status_var.set(f"Ошибка: {exc}")
        messagebox.showerror("Ошибка анализа", str(exc), parent=self.root)

    def _run_background(self, work, on_success):
        """Run data work off Tk's thread; deliver results only from the Tk timer."""
        if self._background_future is not None:
            raise RuntimeError("Расчёт уже выполняется.")
        future = Future()
        self._background_future = future

        def run():
            try:
                future.set_result(work())
            except BaseException as exc:
                future.set_exception(exc)

        def poll():
            self._preload_after_id = None
            if not future.done():
                self._preload_after_id = self.root.after(50, poll)
                return
            self._background_future = None
            try:
                on_success(future.result())
            except Exception:
                self.report_callback_exception(*sys.exc_info())

        Thread(target=run, daemon=True, name="omega-worker").start()
        self._preload_after_id = self.root.after(50, poll)

    def show_batch_progress_window(self):
        self.close_batch_progress_window()
        window = tk.Toplevel(self.root)
        window.title("Анализ батча")
        window.geometry("520x165")
        window.resizable(False, False)
        window.transient(self.root)
        window.protocol("WM_DELETE_WINDOW", lambda: None)

        frame = ttk.Frame(window, padding=20)
        frame.pack(fill="both", expand=True)
        ttk.Label(frame, textvariable=self.batch_progress_label_var, font=("Segoe UI", 12, "bold")).pack(anchor="w")
        ttk.Label(frame, textvariable=self.batch_progress_detail_var).pack(anchor="w", pady=(8, 14))
        self.batch_progress_bar = ttk.Progressbar(frame, mode="determinate", maximum=max(len(self.loaded_batches), 1))
        self.batch_progress_bar.pack(fill="x")
        ttk.Label(frame, text="Дождитесь завершения анализа всех проб.").pack(anchor="w", pady=(12, 0))

        self.batch_progress_window = window
        window.grab_set()
        window.lift()
        window.focus_force()

    def close_batch_progress_window(self):
        window = self.batch_progress_window
        self.batch_progress_window = None
        if self.batch_progress_bar is not None:
            self.batch_progress_bar.stop()
        self.batch_progress_bar = None
        if window is not None and window.winfo_exists():
            try:
                window.grab_release()
            except tk.TclError:
                pass
            window.destroy()

    def preload_loaded_batches(self):
        self._preload_after_id = None
        if not self.loaded_batches:
            return
        total = len(self.loaded_batches)
        while (
            self._preload_batch_index < total
            and self.loaded_batches[self._preload_batch_index].get("processed_df") is not None
        ):
            self._preload_batch_index += 1

        if self._preload_batch_index >= total:
            if self.batch_progress_bar is not None:
                self.batch_progress_bar["value"] = total
            self.close_batch_progress_window()
            self.load_batch_at_index(0)
            self.populate_main_batch_tree()
            if self.batch_results_window is not None and self.batch_results_window.winfo_exists():
                self.populate_batch_results_tree()
            self.status_var.set(f"Рассчёт всех проб завершён: {total}/{total}")
            logging.getLogger("omega").info("CSV analysis and display complete: %d samples", total)
            return

        index = self._preload_batch_index
        batch = self.loaded_batches[index]
        sample_name = batch.get("sample_name", batch.get("file_name", f"Проба {index + 1}"))
        logging.getLogger("omega").info("Processing sample %d/%d: %s", index + 1, total, sample_name)
        self.batch_progress_label_var.set(f"Проба {index + 1} из {total}")
        self.batch_progress_detail_var.set(f"Сейчас анализируется: {sample_name}")
        if self.batch_progress_bar is not None:
            self.batch_progress_bar["maximum"] = total
            self.batch_progress_bar["value"] = index
        self.root.update_idletasks()
        self.status_var.set(f"Расчёт пробы: {index + 1}/{total}")
        reference_targets = self.reference_targets.copy(deep=True)

        def finished(result):
            batch.update(result)
            self._preload_batch_index += 1
            self._preload_after_id = self.root.after(25, self.preload_loaded_batches)

        self._run_background(
            lambda: process_chromatogram_batch(batch["dataframe"], reference_targets),
            finished,
        )

    def refresh_peaks(self, preserve_plot_view: bool = False, redraw_plot: bool = True):
        if self.df_processed is None:
            return
        current_batch = self.loaded_batches[self.current_batch_index] if self.loaded_batches else None
        baseline_mode = current_batch.get("baseline_mode", "chebyshev") if current_batch is not None else "chebyshev"
        scaled_result = finalize_result({**(current_batch or {}),
            "processed_df": self.df_processed, "matched_targets_df": self.matched_targets_df,
            "peaks_df": self.peaks_df, "baseline_mode": baseline_mode}, self.active_instrument_profile)
        self.matched_targets_df = scaled_result["matched_targets_df"]
        cluster_quality_score = scaled_result["cluster_quality_score"]
        confidence = scaled_result["confidence"]
        omega = scaled_result["omega"]
        report_value = scaled_result["omega_report"]
        if current_batch is not None:
            current_batch.update(scaled_result)
            current_batch["best_window"] = self.best_window
            current_batch["rt_shift"] = self.current_rt_shift

        if np.isfinite(report_value):
            qualifier = " (оценка по модели)" if scaled_result["report_status"] == "model_estimate" else ""
            self.omega_var.set(
                f"{report_value:.2f}%{qualifier}"
            )
        else:
            self.omega_var.set("—")

        gamma_text = "γ-Linolenic: —"
        if not self.matched_targets_df.empty:
            gamma_match = self.matched_targets_df[self.matched_targets_df["code"] == "C18:3N6"]
            if not gamma_match.empty:
                gamma_area = float(pd.to_numeric(gamma_match["area"], errors="coerce").iloc[0])
                detected_total_area = float(pd.to_numeric(self.peaks_df["area"], errors="coerce").fillna(0.0).sum())
                gamma_percent = (100.0 * gamma_area / detected_total_area) if detected_total_area > 0 else np.nan
                if np.isfinite(gamma_area):
                    gamma_text = f"γ-Linolenic: area {gamma_area:.2f}"
                    if np.isfinite(gamma_percent):
                        gamma_text = f"{gamma_text} | peaks {gamma_percent:.2f}%"
        self.gamma_var.set(gamma_text)
        self.current_confidence = confidence
        self.confidence_var.set(confidence.get("button_text", "Качество пиков: —"))
        self.confidence_button.state(["!disabled"])
        self.update_review_panel()

        if redraw_plot:
            self.update_plot(preserve_view=preserve_plot_view)
        self.update_table()

        local_shifts = pd.to_numeric(
            self.matched_targets_df.get("rt_local_shift", pd.Series(dtype=float)), errors="coerce"
        ).dropna()
        if not local_shifts.empty and bool(self.matched_targets_df.get(
            "instrument_profile_custom_rt", pd.Series(dtype=bool)
        ).fillna(False).any()):
            rt_text = f"RT correction: {local_shifts.iloc[0]:+.3f}…{local_shifts.iloc[-1]:+.3f} min"
        else:
            rt_text = f"RT shift: {self.current_rt_shift:+.3f} min"
        self.status_var.set(
            f"Определено кислот: {int(self.matched_targets_df['found_rt'].notna().sum())} из {len(self.reference_targets)} · Выберите кислоту для проверки границ"
        )
        self.integration_var.set(
            f"Integration: {len(self.peaks_df)} peaks | SG {self.best_window}"
        )
        self.populate_main_batch_tree()
        if self.batch_results_window is not None and self.batch_results_window.winfo_exists():
            self.populate_batch_results_tree()


    def update_table(self):
        for i in self.tree.get_children():
            self.tree.delete(i)
        if self.matched_targets_df.empty:
            self.selected_target_code = None
            return
        risk = (self.current_confidence or {}).get("high_error_risk", {})
        risky_codes = set(risk.get("peak_codes", [])) if isinstance(risk, dict) else set()
        available_codes = set()
        for _, row in self.matched_targets_df.iterrows():
            code = str(row.get("code", ""))
            available_codes.add(code)
            display_name = str(row.get("display_name", ""))
            if code in risky_codes:
                display_name = f"⚠ {display_name}"
            area = pd.to_numeric(pd.Series([row.get("area")]), errors="coerce").iloc[0]
            area_text = f"{float(area):,.1f}".replace(",", " ") if np.isfinite(area) else "—"
            vals = (
                display_name,
                area_text,
                "" if pd.isna(row.get("percent_area")) else f"{row['percent_area']:.2f}",
                row.get("code", ""),
                "" if pd.isna(row.get("corrected_target_rt", row.get("expected_rt")))
                else f"{float(row.get('corrected_target_rt', row.get('expected_rt'))):.4f}",
                "" if pd.isna(row.get("found_rt")) else f"{row['found_rt']:.4f}",
                readable_status(row.get("status", "")),
            )
            self.tree.insert("", "end", iid=code, values=vals, tags=("review",) if code in risky_codes else ())
        if self.selected_target_code not in available_codes:
            self.selected_target_code = None
        if self.selected_target_code is not None:
            self.tree.selection_set(self.selected_target_code)
            self.tree.focus(self.selected_target_code)


if __name__ == "__main__":
    root = tk.Tk()
    app = ChromatogramApp(root)
    root.after(150, app.start_application)
    root.mainloop()
