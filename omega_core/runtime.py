"""The numerical runtime is part of the analytical method."""
from functools import lru_cache
import hashlib
from importlib import metadata
import json
import platform
from pathlib import Path
import sys

VERSIONS = {"numpy": "2.5.3", "pandas": "2.3.3", "scipy": "1.18.1",
            "pybaselines": "1.2.1", "lmfit": "1.3.4", "pyopenms": "3.5.0"}


def installed_versions():
    versions = {}
    for name in VERSIONS:
        try:
            versions[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            versions[name] = None
    return versions


@lru_cache(maxsize=1)
def verify_versions():
    installed = installed_versions()
    differences = [f"{name}: {installed[name] or 'не установлен'}; требуется {expected}"
                   for name, expected in VERSIONS.items() if installed[name] != expected]
    if differences:
        raise RuntimeError("Окружение расчёта отличается от проверенного.\n" + "\n".join(differences)
                           + "\nИспользуйте одну portable-сборку или установите requirements.txt в окружение проекта.")


@lru_cache(maxsize=1)
def code_digest():
    root = Path(__file__).resolve().parent
    manifest = root.parent / "runtime_manifest.json"
    if getattr(sys, "frozen", False) and manifest.exists():
        return json.loads(manifest.read_text(encoding="utf-8"))["engine_sha256"]
    digest = hashlib.sha256()
    for path in sorted(root.glob("*.py")):
        digest.update(path.name.encode())
        digest.update(path.read_text(encoding="utf-8").encode())
    return digest.hexdigest()


def snapshot(reference_targets=None):
    from . import signal, clusters, metrics, chromatopy_adapter, rt_profile, joint_clusters, peak_geometry, assignment
    settings = {}
    for module in (signal, clusters, metrics, chromatopy_adapter, rt_profile, joint_clusters, peak_geometry, assignment):
        for name, value in vars(module).items():
            if name.startswith("_") or not name.isupper():
                continue
            if hasattr(value, "tolist"):
                value = value.tolist()
            elif isinstance(value, (set, frozenset)):
                value = sorted(value)
            try:
                json.dumps(value, allow_nan=False)
            except (ValueError, TypeError):
                continue
            settings[f"{module.__name__}.{name}"] = value
    report = {"python": platform.python_version(), "platform": platform.platform(),
              "versions": installed_versions(), "engine_sha256": code_digest(), "settings": settings}
    if reference_targets is not None:
        report["targets_sha256"] = hashlib.sha256(reference_targets.to_json(orient="split", double_precision=15).encode()).hexdigest()
    return report
