"""Manual peak selection, boundary editing and edit import/export."""
import tkinter as tk
from tkinter import filedialog, messagebox
import numpy as np
import pandas as pd
from omega_core import integration, manual_edits
from omega_core.results import finalize_result
from omega_plot_ui import _get_x_column_name


class ManualIntegrationUI:
    def handle_target_selection(self, event=None):
        self.cancel_manual_pick()
        selection = self.tree.selection()
        code = selection[0] if selection else None
        if code == self.selected_target_code:
            return
        self.selected_target_code = code
        self.load_selected_integration_bounds(silent=True)
        self.update_peak_details()
        if self.df_processed is not None:
            group = next((g for g in ('C16','C18','C20','C22') if str(code).startswith(g+':')), 'Пик')
            self.view_var.set(group)
            self.change_plot_view()

    def load_selected_integration_bounds(self, silent: bool = False):
        if not self.selected_target_code or self.matched_targets_df.empty:
            self.manual_start_var.set("")
            self.manual_end_var.set("")
            if not silent:
                self.status_var.set("Выбери пик в таблице перед ручной интеграцией.")
            return

        row = self.matched_targets_df[self.matched_targets_df["code"] == self.selected_target_code]
        if row.empty:
            self.manual_start_var.set("")
            self.manual_end_var.set("")
            return
        start_x = pd.to_numeric(row["integration_start_x"], errors="coerce").iloc[0]
        end_x = pd.to_numeric(row["integration_end_x"], errors="coerce").iloc[0]
        self.manual_start_var.set("" if not np.isfinite(start_x) else f"{float(start_x):.8f}")
        self.manual_end_var.set("" if not np.isfinite(end_x) else f"{float(end_x):.8f}")
        if not silent:
            self.status_var.set(f"Границы {self.selected_target_code} загружены для ручной правки.")

    def apply_manual_integration(self):
        if self.df_processed is None or self.matched_targets_df.empty:
            messagebox.showwarning("Ручная интеграция", "Сначала открой CSV и выбери образец.", parent=self.root)
            return
        if not self.selected_target_code:
            messagebox.showwarning("Ручная интеграция", "Сначала выбери пик в таблице справа.", parent=self.root)
            return

        try:
            start_x = float(str(self.manual_start_var.get()).replace(",", "."))
            end_x = float(str(self.manual_end_var.get()).replace(",", "."))
        except ValueError:
            messagebox.showerror("Ручная интеграция", "Начало и конец интервала должны быть числами.", parent=self.root)
            return
        if not (np.isfinite(start_x) and np.isfinite(end_x) and end_x > start_x):
            messagebox.showerror("Ручная интеграция", "Конец интервала должен быть больше начала.", parent=self.root)
            return

        batch = self.loaded_batches[self.current_batch_index]
        try:
            result = manual_edits.apply_edit(batch, self.selected_target_code, start_x, end_x)
        except ValueError as exc:
            messagebox.showerror("Ручная интеграция", str(exc), parent=self.root)
            return
        self._refresh_after_manual_edit()
        self.status_var.set(f"Ручная интеграция {self.selected_target_code}: {result.start:.5f}–{result.end:.5f}, площадь {result.area:.4f}")

    def _refresh_after_manual_edit(self):
        batch = self.loaded_batches[self.current_batch_index]
        self.matched_targets_df = batch["matched_targets_df"]
        self.load_selected_integration_bounds(silent=True)
        # Redraw from the committed result, including its exact apex. The fast
        # drag overlay is only a preview and must not survive apply/undo.
        self.refresh_peaks(preserve_plot_view=True)

    def propose_signal_bounds(self):
        if self.df_processed is None or not self.selected_target_code:
            self.status_var.set("Выберите пик в таблице.")
            return
        row = self.matched_targets_df[self.matched_targets_df["code"] == self.selected_target_code].iloc[0]
        try:
            from omega_core.peak_geometry import proposal_for_apex
            proposal = proposal_for_apex(self.peaks_df, float(row['found_rt']))
        except ValueError as exc:
            self.status_var.set(str(exc))
            return
        if not proposal.resolved:
            self.status_var.set("Границы неоднозначны: нет устойчивой впадины или выхода на базовую линию. Нужна ручная проверка.")
            return
        self.manual_start_var.set(f"{proposal.start:.8f}")
        self.manual_end_var.set(f"{proposal.end:.8f}")
        self.update_plot(preserve_view=True)
        names = {"baseline": "базовая линия", "valley": "впадина"}
        self.status_var.set(f"Предложение: слева {names[proposal.left_kind]}, справа {names[proposal.right_kind]}. Проверьте и нажмите «Применить».")

    def undo_manual_edit(self):
        if self.loaded_batches and manual_edits.undo_edit(self.loaded_batches[self.current_batch_index]):
            self._refresh_after_manual_edit()
            self.status_var.set("Последняя правка этой пробы отменена.")

    def save_manual_edits(self):
        if not self.current_file or not self.loaded_batches:
            return
        path = filedialog.asksaveasfilename(parent=self.root, defaultextension=".json", initialfile="omega_manual_edits.json")
        if path:
            manual_edits.export_edits(path, self.current_file, self.active_instrument_profile, self.loaded_batches)
            self.status_var.set("Ручные правки сохранены.")
            return True
        return False

    def load_manual_edits(self):
        if not self.current_file or not self.loaded_batches:
            return
        path = filedialog.askopenfilename(parent=self.root, filetypes=[("Правки Omega", "*.json")])
        if path:
            try:
                count = manual_edits.import_edits(path, self.current_file, self.active_instrument_profile, self.loaded_batches)
            except (ValueError, KeyError, TypeError, OSError) as exc:
                messagebox.showerror("Загрузка правок", str(exc), parent=self.root)
                return
            # Refresh every changed result, not just the displayed sample.
            for batch in self.loaded_batches:
                if batch.get("manual_overrides"):
                    batch.update(finalize_result(batch, self.active_instrument_profile))
            self._refresh_after_manual_edit()
            self.populate_main_batch_tree()
            self.status_var.set(f"Восстановлены правки: проб {count}.")

    def _manual_bounds_from_vars(self):
        try:
            start_x = float(str(self.manual_start_var.get()).replace(",", "."))
            end_x = float(str(self.manual_end_var.get()).replace(",", "."))
        except ValueError:
            return np.nan, np.nan
        return start_x, end_x

    def _selected_manual_drag_bounds(self):
        if not self.selected_target_code or self.matched_targets_df.empty:
            return np.nan, np.nan
        start_x, end_x = self._manual_bounds_from_vars()
        if np.isfinite(start_x) and np.isfinite(end_x):
            return start_x, end_x
        row = self.matched_targets_df[self.matched_targets_df["code"] == self.selected_target_code]
        if row.empty:
            return np.nan, np.nan
        start_x = pd.to_numeric(row["integration_start_x"], errors="coerce").iloc[0]
        end_x = pd.to_numeric(row["integration_end_x"], errors="coerce").iloc[0]
        return float(start_x), float(end_x)

    def cancel_manual_pick(self):
        self._insert_pick = None
        self._manual_pick = None
        if self._manual_pick_artist is not None:
            try:
                self._manual_pick_artist.remove()
            except (ValueError, NotImplementedError):
                pass
            self._manual_pick_artist = None
            self.canvas.draw_idle()
        if hasattr(self, "manual_pick_button"):
            self.manual_pick_button.configure(text="Выделить границы")

    def toggle_manual_pick(self):
        if self._manual_pick is not None:
            self.cancel_manual_pick()
            self.status_var.set("Установка границ отменена.")
            return
        if self.df_processed is None or not self.selected_target_code:
            self.status_var.set("Сначала выберите компонент в таблице.")
            return
        if self.plot_toolbar.mode == "pan/zoom":
            self.plot_toolbar.pan()
        elif self.plot_toolbar.mode == "zoom rect":
            self.plot_toolbar.zoom()
        self._manual_pick = {"code": self.selected_target_code, "batch": self.current_batch_index, "first": None}
        self.manual_pick_button.configure(text="Отменить выбор (Esc)")
        self.status_var.set("Щёлкните две границы пика на графике. Второй щелчок применит интервал; Esc отменяет.")

    def handle_manual_boundary_press(self, event):
        if self.insertion_click(event):
            return
        if event.button != 1 or event.inaxes not in [self.ax, *getattr(self, "preview_axes", [])]:
            return
        if event.xdata is None or self.df_processed is None or not self.selected_target_code:
            return
        if self._manual_pick is None and self.plot_toolbar.mode:
            return
        if self._manual_pick is not None:
            pick = self._manual_pick
            if pick["code"] != self.selected_target_code or pick["batch"] != self.current_batch_index:
                self.cancel_manual_pick()
                return
            x = self.df_processed[_get_x_column_name(self.df_processed)].to_numpy()
            value = float(event.xdata)
            if not np.isfinite(value) or not x[0] <= value <= x[-1]:
                return
            if pick["first"] is None:
                pick["first"] = value
                self._manual_pick_artist = event.inaxes.axvline(value, color="#2463C4", linestyle="--")
                self.canvas.draw_idle()
                self.status_var.set("Теперь щёлкните вторую границу пика.")
            elif value != pick["first"]:
                start, end = sorted((value, pick["first"]))
                self.cancel_manual_pick()
                self.manual_start_var.set(f"{start:.8f}")
                self.manual_end_var.set(f"{end:.8f}")
                self.apply_manual_integration()
            return
        start_x, end_x = self._selected_manual_drag_bounds()
        if not (np.isfinite(start_x) and np.isfinite(end_x)):
            return

        x_min, x_max = event.inaxes.get_xlim()
        tolerance = max(0.0025, abs(float(x_max) - float(x_min)) * 0.010)
        distances = {"start": abs(float(event.xdata) - start_x), "end": abs(float(event.xdata) - end_x)}
        boundary, distance = min(distances.items(), key=lambda item: item[1])
        if distance > tolerance:
            return
        self._manual_drag_active_boundary = boundary
        self._manual_drag_axis = event.inaxes
        self.status_var.set(f"Тяни {'левую' if boundary == 'start' else 'правую'} границу {self.selected_target_code} мышкой…")

    def handle_manual_boundary_motion(self, event):
        if self._manual_drag_active_boundary is None or event.xdata is None or self.df_processed is None:
            return
        x_col = _get_x_column_name(self.df_processed)
        x_values = self.df_processed[x_col].to_numpy(dtype=float)
        x_value = float(np.clip(float(event.xdata), float(np.nanmin(x_values)), float(np.nanmax(x_values))))
        start_x, end_x = self._selected_manual_drag_bounds()
        if self._manual_drag_active_boundary == "start":
            if np.isfinite(end_x):
                x_value = min(x_value, float(end_x) - 1e-5)
            self.manual_start_var.set(f"{x_value:.8f}")
        else:
            if np.isfinite(start_x):
                x_value = max(x_value, float(start_x) + 1e-5)
            self.manual_end_var.set(f"{x_value:.8f}")
        self._manual_drag_pending_bounds = self._manual_bounds_from_vars()
        if self._manual_drag_after_id is None:
            self._manual_drag_after_id = self.root.after(16, self._flush_manual_drag_overlay)

    def _flush_manual_drag_overlay(self):
        self._manual_drag_after_id = None
        bounds = self._manual_drag_pending_bounds
        self._manual_drag_pending_bounds = None
        if bounds is None:
            return
        start_x, end_x = bounds
        if not np.all(np.isfinite([start_x, end_x])) or end_x <= start_x:
            return
        drag_axis = self._manual_drag_axis
        overlay_items = (
            [(drag_axis, self._manual_overlay_artists[drag_axis])]
            if drag_axis in self._manual_overlay_artists
            else list(self._manual_overlay_artists.items())
        )
        changed_axes = []
        for axis, artists in overlay_items:
            try:
                x = artists["x"]
                fill_y = artists["fill_y"]
                start_idx = max(0, min(int(np.searchsorted(x, start_x, side="left")), len(x) - 2))
                end_idx = max(start_idx + 1, min(int(np.searchsorted(x, end_x, side="right") - 1), len(x) - 1))
                artists["fill"].remove()
                fill_x, exact_y = integration.interval_points(x, fill_y, start_x, end_x)
                artists["fill"] = axis.fill_between(
                    fill_x,
                    0.0,
                    np.maximum(exact_y, 0.0),
                    color="#2463C4",
                    alpha=0.36 if artists["compact"] else 0.40,
                    linewidth=0.0,
                    zorder=3,
                )
                artists["start_line"].set_xdata([start_x, start_x])
                artists["end_line"].set_xdata([end_x, end_x])
                apex_idx = int(start_idx + np.argmax(fill_y[start_idx:end_idx + 1]))
                marker = artists.get("marker")
                if marker is not None:
                    marker.set_offsets(np.asarray([[x[apex_idx], artists["marker_y"][apex_idx]]]))
                annotation = artists.get("annotation")
                if annotation is not None:
                    annotation.xy = (x[apex_idx], artists["marker_y"][apex_idx])
                    annotation.set_text(f"{self.selected_target_code}  RT {x[apex_idx]:.4f}")
                changed_axes.append(axis)
            except (KeyError, ValueError, RuntimeError):
                continue
        try:
            renderer = self.canvas.get_renderer()
            for axis in changed_axes:
                axis.draw(renderer)
                self.canvas.blit(axis.bbox)
        except (AttributeError, RuntimeError, tk.TclError):
            self.canvas.draw_idle()

    def handle_manual_boundary_release(self, event):
        if self._manual_drag_active_boundary is None:
            return
        self._manual_drag_active_boundary = None
        if self._manual_drag_after_id is not None:
            self.root.after_cancel(self._manual_drag_after_id)
            self._manual_drag_after_id = None
        self._manual_drag_pending_bounds = None
        self._manual_drag_axis = None
        self.status_var.set(f"Граница {self.selected_target_code} изменена мышью; пересчитываю пик…")
        self.apply_manual_integration()
