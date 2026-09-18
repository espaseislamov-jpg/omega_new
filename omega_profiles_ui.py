"""Instrument profile selection, editing and recalculation controls."""
import json
from pathlib import Path
import tkinter as tk
from tkinter import filedialog, messagebox, simpledialog, ttk
import numpy as np
import pandas as pd
from omega_core import instrument_profiles


class ProfileManagerUI:
    def _refresh_profile_combo(self):
        if self.profile_combo is None:
            return
        names = [profile["name"] for profile in self.profile_store.get("profiles", [])]
        self.profile_combo.configure(values=names)
        self.profile_var.set(self.active_instrument_profile["name"])

    def _profile_by_id(self, profile_id: str):
        for profile in self.profile_store.get("profiles", []):
            if str(profile.get("id")) == str(profile_id):
                return profile
        return None

    def on_profile_combo_selected(self, event=None):
        name = self.profile_var.get()
        profile = next(
            (item for item in self.profile_store.get("profiles", []) if item.get("name") == name),
            None,
        )
        if profile is None:
            self.profile_var.set(self.active_instrument_profile["name"])
            return
        self.activate_instrument_profile(profile["id"])

    def activate_instrument_profile(self, profile_id: str, force_recalculate: bool = False):
        profile = self._profile_by_id(profile_id)
        if profile is None:
            return False
        if self._new_profile_dialog is not None and self._new_profile_dialog.winfo_exists():
            self.new_profile_var.set(profile['name'])
            return True
        if self._background_future is not None or self._preload_after_id is not None:
            self.profile_var.set(self.active_instrument_profile['name'])
            return False
        changed = str(profile_id) != str(self.active_instrument_profile.get("id"))
        if self.loaded_batches and (changed or force_recalculate):
            confirmed = messagebox.askyesno(
                "Сменить профиль прибора",
                "Все загруженные пробы будут последовательно пересчитаны с новым профилем. Продолжить?",
                parent=self.profile_window if self.profile_window is not None else self.root,
            )
            if not confirmed:
                self.profile_var.set(self.active_instrument_profile["name"])
                return False

        self.profile_store["active_profile_id"] = str(profile_id)
        self.profile_store = instrument_profiles.save_store(self.profile_store, self.base_reference_targets)
        self.active_instrument_profile = instrument_profiles.active_profile(self.profile_store)
        self.reference_targets = instrument_profiles.apply_profile_to_targets(
            self.base_reference_targets, self.active_instrument_profile
        )
        self._refresh_profile_combo()
        if self.loaded_batches and (changed or force_recalculate):
            self.recalculate_loaded_batches()
        else:
            self.status_var.set(f"Активный профиль: {self.active_instrument_profile['name']}")
        return True

    def recalculate_loaded_batches(self):
        if self._preload_after_id is not None:
            self.root.after_cancel(self._preload_after_id)
            self._preload_after_id = None
        preserved = {"file_name", "signal_name", "acquired_at", "source_path", "sample_name", "dataframe"}
        for batch in self.loaded_batches:
            for key in list(batch):
                if key not in preserved:
                    del batch[key]
        self.df_processed = None
        self.peaks_df = pd.DataFrame()
        self.matched_targets_df = pd.DataFrame()
        self._preload_batch_index = 0
        self.show_batch_progress_window()
        self._preload_after_id = self.root.after(50, self.preload_loaded_batches)

    def open_profile_manager(self):
        if self.profile_window is not None and self.profile_window.winfo_exists():
            self.profile_window.deiconify()
            self.profile_window.lift()
            return

        window = tk.Toplevel(self.root)
        self.profile_window = window
        window.title("Профили приборов и колонок")
        window.geometry("930x620")
        window.minsize(760, 500)
        window.transient(self.root)

        container = ttk.Frame(window, padding=12)
        container.pack(fill="both", expand=True)
        panes = ttk.Panedwindow(container, orient="horizontal")
        panes.pack(fill="both", expand=True)
        left = ttk.Frame(panes)
        right = ttk.Frame(panes, padding=(12, 0, 0, 0))
        panes.add(left, weight=2)
        panes.add(right, weight=3)

        profile_tree = ttk.Treeview(
            left,
            columns=("name", "multiplier", "judge"),
            show="headings",
            selectmode="browse",
        )
        profile_tree.heading("name", text="Профиль")
        profile_tree.heading("multiplier", text="Множитель")
        profile_tree.heading("judge", text="Судья")
        profile_tree.column("name", width=150, anchor="w")
        profile_tree.column("multiplier", width=65, anchor="center", stretch=False)
        profile_tree.column("judge", width=165, anchor="w")
        profile_tree.pack(fill="both", expand=True)

        left_buttons = ttk.Frame(left)
        left_buttons.pack(fill="x", pady=(8, 0))

        profile_name_var = tk.StringVar()
        multiplier_var = tk.StringVar(value="1.0")
        calculation_labels = {
            "Историческая формула": instrument_profiles.LEGACY_CALCULATION_MODE,
            "Прямой расчёт по выбранным площадям": instrument_profiles.DIRECT_COMPONENTS_CALCULATION_MODE,
        }
        calculation_mode_var = tk.StringVar(value="Прямой расчёт по выбранным площадям")
        component_codes_var = tk.StringVar(
            value=", ".join(instrument_profiles.DEFAULT_OMEGA_COMPONENT_CODES)
        )
        judge_status_var = tk.StringVar(value="Судья: —")
        area_baseline_labels = {"По геометрии пиков": "geometry", "Chebyshev": "chebyshev"}
        area_baseline_var = tk.StringVar(value="Chebyshev")
        form = ttk.Frame(right)
        form.pack(fill="x")
        ttk.Label(form, text="Название").grid(row=0, column=0, sticky="w")
        name_entry = ttk.Entry(form, textvariable=profile_name_var)
        name_entry.grid(row=0, column=1, sticky="ew", padx=(8, 0))
        ttk.Label(form, text="Множитель итогового результата").grid(row=1, column=0, sticky="w", pady=(8, 0))
        multiplier_entry = ttk.Entry(form, textvariable=multiplier_var, width=14)
        multiplier_entry.grid(row=1, column=1, sticky="w", padx=(8, 0), pady=(8, 0))
        ttk.Label(form, text="Схема расчёта").grid(row=2, column=0, sticky="w", pady=(8, 0))
        calculation_combo = ttk.Combobox(
            form,
            textvariable=calculation_mode_var,
            values=list(calculation_labels),
            state="readonly",
            width=38,
        )
        calculation_combo.grid(row=2, column=1, sticky="ew", padx=(8, 0), pady=(8, 0))
        ttk.Label(form, text="Пики числителя через запятую").grid(
            row=3, column=0, sticky="w", pady=(8, 0)
        )
        component_codes_entry = ttk.Entry(form, textvariable=component_codes_var)
        component_codes_entry.grid(row=3, column=1, sticky="ew", padx=(8, 0), pady=(8, 0))
        ttk.Label(form, text="Калибровка судьи").grid(row=4, column=0, sticky="w", pady=(8, 0))
        ttk.Label(form, textvariable=judge_status_var).grid(
            row=4, column=1, sticky="w", padx=(8, 0), pady=(8, 0)
        )
        form.columnconfigure(1, weight=1)
        ttk.Label(form, text="Базовая линия площадей").grid(row=5, column=0, sticky="w", pady=(8, 0))
        area_baseline_combo = ttk.Combobox(form, textvariable=area_baseline_var,
                                         values=list(area_baseline_labels), state="readonly")
        area_baseline_combo.grid(row=5, column=1, sticky="ew", padx=(8, 0), pady=(8, 0))

        ttk.Label(
            right,
            text="Времена выхода, мин. Двойной щелчок по строке — изменить RT.",
        ).pack(fill="x", pady=(12, 5))
        rt_tree = ttk.Treeview(right, columns=("code", "display", "rt"), show="headings", selectmode="browse")
        rt_tree.heading("code", text="Код")
        rt_tree.heading("display", text="Кислота")
        rt_tree.heading("rt", text="RT, мин")
        rt_tree.column("code", width=100, anchor="center", stretch=False)
        rt_tree.column("display", width=245, anchor="w")
        rt_tree.column("rt", width=100, anchor="center", stretch=False)
        rt_tree.pack(fill="both", expand=True)

        actions = ttk.Frame(right)
        actions.pack(fill="x", pady=(8, 0))

        def selected_id():
            selection = profile_tree.selection()
            return selection[0] if selection else None

        def load_details(event=None):
            profile = self._profile_by_id(selected_id()) if selected_id() else None
            if profile is None:
                return
            profile_name_var.set(profile["name"])
            multiplier_var.set(f"{float(profile.get('result_multiplier', 1.0)):.8g}")
            mode = str(
                profile.get(
                    "calculation_mode", instrument_profiles.LEGACY_CALCULATION_MODE
                )
            )
            calculation_mode_var.set(
                next(
                    (label for label, value in calculation_labels.items() if value == mode),
                    "Историческая формула",
                )
            )
            component_codes_var.set(
                ", ".join(
                    profile.get(
                        "omega_component_codes",
                        instrument_profiles.DEFAULT_OMEGA_COMPONENT_CODES,
                    )
                )
            )
            judge_status_var.set(instrument_profiles.judge_calibration_label(profile))
            area_baseline_var.set(next(label for label,value in area_baseline_labels.items()
                                       if value==profile.get('integration_baseline','geometry')))
            for item in rt_tree.get_children():
                rt_tree.delete(item)
            names = self.base_reference_targets.set_index("code")["display_name"].to_dict()
            for code in self.base_reference_targets.sort_values("order_index")["code"].astype(str):
                value = profile["retention_times"].get(code, np.nan)
                rt_tree.insert("", "end", iid=code, values=(code, names.get(code, code), f"{float(value):.4f}"))
            locked = str(profile.get("id")) == instrument_profiles.LEGACY_PROFILE_ID
            state = "disabled" if locked else "normal"
            name_entry.configure(state=state)
            multiplier_entry.configure(state=state)
            calculation_combo.configure(state="disabled" if locked else "readonly")
            area_baseline_combo.configure(state="disabled" if locked else "readonly")
            component_codes_entry.configure(state=state)

        def reload_profiles(select_profile_id=None):
            for item in profile_tree.get_children():
                profile_tree.delete(item)
            active_id = str(self.profile_store.get("active_profile_id"))
            for profile in self.profile_store.get("profiles", []):
                marker = "● " if str(profile["id"]) == active_id else ""
                profile_tree.insert(
                    "", "end", iid=str(profile["id"]),
                    values=(
                        marker + profile["name"],
                        f"{float(profile.get('result_multiplier', 1.0)):.6g}",
                        instrument_profiles.judge_calibration_label(profile),
                    ),
                )
            target = str(select_profile_id or active_id)
            if profile_tree.exists(target):
                profile_tree.selection_set(target)
                profile_tree.focus(target)
                profile_tree.see(target)
            load_details()

        def edit_rt(event=None):
            profile = self._profile_by_id(selected_id()) if selected_id() else None
            selection = rt_tree.selection()
            if profile is None or not selection:
                return
            if str(profile.get("id")) == instrument_profiles.LEGACY_PROFILE_ID:
                messagebox.showinfo("Встроенный профиль", "Сначала создайте копию встроенного профиля.", parent=window)
                return
            code = selection[0]
            current = float(rt_tree.item(code, "values")[2])
            value = simpledialog.askfloat(
                "Время выхода", f"RT для {code}, мин:", initialvalue=current,
                minvalue=0.001, maxvalue=100.0, parent=window,
            )
            if value is not None:
                values = list(rt_tree.item(code, "values"))
                values[2] = f"{value:.4f}"
                rt_tree.item(code, values=values)

        def save_edits():
            profile_id = selected_id()
            profile = self._profile_by_id(profile_id) if profile_id else None
            if profile is None:
                return False
            if str(profile_id) == instrument_profiles.LEGACY_PROFILE_ID:
                messagebox.showinfo("Встроенный профиль", "Встроенный профиль неизменяем. Создайте его копию.", parent=window)
                return False
            candidate = dict(profile)
            candidate["name"] = profile_name_var.get().strip()
            candidate["result_multiplier"] = multiplier_var.get().strip()
            candidate["calculation_mode"] = calculation_labels.get(
                calculation_mode_var.get(), instrument_profiles.LEGACY_CALCULATION_MODE
            )
            candidate['integration_baseline'] = area_baseline_labels[area_baseline_var.get()]
            candidate["omega_component_codes"] = [
                code.strip()
                for code in component_codes_var.get().split(",")
                if code.strip()
            ]
            candidate["retention_times"] = {
                code: rt_tree.item(code, "values")[2] for code in rt_tree.get_children()
            }
            # Preserve full precision when the operator did not edit a cell.
            for code, value in candidate["retention_times"].items():
                original = float(profile["retention_times"][code])
                if float(value) == float(f"{original:.4f}"):
                    candidate["retention_times"][code] = original
            if any(float(value) != float(profile["retention_times"][code])
                   for code, value in candidate["retention_times"].items()):
                candidate["custom_rt"] = True
            try:
                candidate = instrument_profiles.validate_profile(
                    candidate, self.base_reference_targets.sort_values("order_index")["code"].astype(str)
                )
                duplicate = next(
                    (item for item in self.profile_store["profiles"]
                     if item["id"] != profile_id and item["name"].casefold() == candidate["name"].casefold()),
                    None,
                )
                if duplicate is not None:
                    raise ValueError("Профиль с таким названием уже существует.")
            except (ValueError, TypeError) as exc:
                messagebox.showerror("Профиль не сохранён", str(exc), parent=window)
                return False
            calibration_inputs_changed = bool(
                float(candidate.get("result_multiplier", 1.0))
                != float(profile.get("result_multiplier", 1.0))
                or candidate.get("calculation_mode") != profile.get("calculation_mode")
                or candidate.get('integration_baseline','geometry') != profile.get('integration_baseline','geometry')
                or candidate.get("omega_component_codes") != profile.get("omega_component_codes")
                or candidate.get("retention_times") != profile.get("retention_times")
            )
            if calibration_inputs_changed:
                candidate.update({
                    "judge_calibrated": False,
                    "judge_emergency_enabled": False,
                    "judge_manual_samples": 0,
                    "judge_error_samples": 0,
                    "judge_validation_batches": 0,
                    "judge_c22_height_priority_threshold": None,
                })
            old_profile = dict(profile)
            old_profile["retention_times"] = dict(profile["retention_times"])
            changed_index = None
            for index, item in enumerate(self.profile_store["profiles"]):
                if str(item["id"]) == str(profile_id):
                    self.profile_store["profiles"][index] = candidate
                    changed_index = index
                    break
            self.profile_store = instrument_profiles.save_store(self.profile_store, self.base_reference_targets)
            if str(profile_id) == str(self.active_instrument_profile.get("id")):
                if not self.activate_instrument_profile(profile_id, force_recalculate=True):
                    if changed_index is not None:
                        self.profile_store["profiles"][changed_index] = old_profile
                        self.profile_store = instrument_profiles.save_store(
                            self.profile_store, self.base_reference_targets
                        )
                    return False
            self._refresh_profile_combo()
            reload_profiles(profile_id)
            return True

        def create_profile(copy_selected=False):
            source = self._profile_by_id(selected_id()) if copy_selected and selected_id() else None
            default_name = f"{source['name']} — копия" if source else "Новый прибор"
            name = simpledialog.askstring("Новый профиль", "Название профиля:", initialvalue=default_name, parent=window)
            if not name:
                return
            name = instrument_profiles.unique_profile_name(self.profile_store, name)
            profile = (
                instrument_profiles.copy_profile(source, name)
                if source else instrument_profiles.new_profile(name, self.base_reference_targets)
            )
            self.profile_store["profiles"].append(profile)
            self.profile_store = instrument_profiles.save_store(self.profile_store, self.base_reference_targets)
            reload_profiles(profile["id"])
            self._refresh_profile_combo()

        def delete_profile():
            profile_id = selected_id()
            if not profile_id or profile_id == instrument_profiles.LEGACY_PROFILE_ID:
                return
            if any(p['id'] == profile_id for p in instrument_profiles.builtin_profiles(self.base_reference_targets)):
                messagebox.showinfo("Базовый профиль", "Базовый профиль нельзя удалить. Можно изменить его настройки или создать копию.", parent=window)
                return
            if str(profile_id) == str(self.profile_store.get("active_profile_id")):
                messagebox.showerror("Нельзя удалить", "Сначала сделайте активным другой профиль.", parent=window)
                return
            profile = self._profile_by_id(profile_id)
            if not messagebox.askyesno("Удалить профиль", f"Удалить «{profile['name']}»?", parent=window):
                return
            self.profile_store["profiles"] = [item for item in self.profile_store["profiles"] if item["id"] != profile_id]
            self.profile_store = instrument_profiles.save_store(self.profile_store, self.base_reference_targets)
            reload_profiles()
            self._refresh_profile_combo()

        def activate_selected():
            profile_id = selected_id()
            if profile_id and self.activate_instrument_profile(profile_id):
                reload_profiles(profile_id)

        def export_selected():
            profile = self._profile_by_id(selected_id()) if selected_id() else None
            if profile is None:
                return
            path = filedialog.asksaveasfilename(
                title="Экспорт профиля", defaultextension=instrument_profiles.PROFILE_FILE_SUFFIX,
                filetypes=[("Профиль Omega", "*.omega-profile.json"), ("JSON", "*.json")],
                initialfile=f"{profile['name']}{instrument_profiles.PROFILE_FILE_SUFFIX}", parent=window,
            )
            if path:
                instrument_profiles.export_profile(
                    profile, Path(path), self.base_reference_targets.sort_values("order_index")["code"].astype(str)
                )

        def import_one():
            path = filedialog.askopenfilename(
                title="Импорт профиля", filetypes=[("Профиль Omega", "*.omega-profile.json *.json"), ("Все файлы", "*.*")],
                parent=window,
            )
            if not path:
                return
            try:
                profile = instrument_profiles.import_profile(Path(path), self.profile_store, self.base_reference_targets)
                self.profile_store["profiles"].append(profile)
                self.profile_store = instrument_profiles.save_store(self.profile_store, self.base_reference_targets)
            except (OSError, ValueError, TypeError, json.JSONDecodeError) as exc:
                messagebox.showerror("Импорт не выполнен", str(exc), parent=window)
                return
            reload_profiles(profile["id"])
            self._refresh_profile_combo()

        profile_tree.bind("<<TreeviewSelect>>", load_details)
        rt_tree.bind("<Double-1>", edit_rt)
        ttk.Button(left_buttons, text="Новый", command=lambda: create_profile(False)).pack(side="left")
        ttk.Button(left_buttons, text="Копия", command=lambda: create_profile(True)).pack(side="left", padx=(6, 0))
        ttk.Button(left_buttons, text="Удалить", command=delete_profile).pack(side="left", padx=(6, 0))
        ttk.Button(actions, text="Изменить RT", command=edit_rt).pack(side="left")
        ttk.Button(actions, text="Сохранить", command=save_edits).pack(side="left", padx=(6, 0))
        ttk.Button(actions, text="Сделать активным", command=activate_selected).pack(side="left", padx=(6, 0))
        ttk.Button(actions, text="Импорт", command=import_one).pack(side="right")
        ttk.Button(actions, text="Экспорт", command=export_selected).pack(side="right", padx=(0, 6))

        def on_close():
            window.destroy()
            self.profile_window = None

        window.protocol("WM_DELETE_WINDOW", on_close)
        reload_profiles()
