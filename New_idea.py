"""Tkinter GUI for the production Omega chromatogram engine."""

import re
from pathlib import Path

from omega_path_compat import configure_windows_path_compat

configure_windows_path_compat()

import numpy as np
import pandas as pd
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

import omega_core
from omega_core import metrics as core_metrics
from omega_core.io import ensure_runtime_file


PREVIEW_WINDOWS = [
    ("6.0-7.3", 6.0, 7.3),
    ("7.4-7.7", 7.4, 7.7),
    ("8.3-8.7", 8.3, 8.7),
    ("9.1-9.4", 9.1, 9.4),
]


def _get_x_column_name(df: pd.DataFrame) -> str:
    if "x_corrected" in df.columns:
        return "x_corrected"
    if "x" in df.columns:
        return "x"
    raise KeyError("Expected an x or x_corrected column")


def _recompute_matched_percent_area(matched_targets_df: pd.DataFrame) -> pd.DataFrame:
    out = matched_targets_df.copy()
    areas = pd.to_numeric(out["area"], errors="coerce")
    total_area = float(areas.fillna(0.0).sum())
    out["percent_area"] = 100.0 * areas / total_area if total_area > 0 else np.nan
    return out


def process_chromatogram_batch(dataframe: pd.DataFrame, reference_targets: pd.DataFrame) -> dict:
    """Run the same modular engine used by the regression harness."""
    engine_input = dataframe.copy()
    if "x_corrected" not in engine_input.columns and "x" in engine_input.columns:
        engine_input["x_corrected"] = pd.to_numeric(engine_input["x"], errors="coerce")
    result = dict(omega_core.process_batch(engine_input, reference_targets))
    result["engine"] = "omega_core"
    result.setdefault("total_area", result.get("omega", {}).get("total_area", np.nan))
    return result


class ChromatogramApp:
    def __init__(self, root: tk.Tk):
        self.root = root
        self.root.title("Chromatogram Peak Detector (Integrated)")
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
        self.reference_targets = omega_core.load_reference_targets(self.reference_json_path)

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
        self._manual_overlay_artists = {}

        self.status_var = tk.StringVar(value="Выбери CSV-файл.")
        self.file_var = tk.StringVar(value="Файл не выбран")
        self.omega_var = tk.StringVar(value="Omega-3: —")
        self.integration_var = tk.StringVar(value="Integration: —")
        self.gamma_var = tk.StringVar(value="γ-Linolenic: —")
        self.batch_var = tk.StringVar(value="Series: —")
        self.confidence_var = tk.StringVar(value="Качество пиков: —")
        self.current_confidence = None

        self._build_ui()

    def _build_ui(self):
        controls = ttk.Frame(self.root, padding=(10, 8))
        controls.pack(fill="x")
        action_bar = ttk.Frame(controls)
        action_bar.pack(fill="x")
        self.confidence_button = ttk.Button(
            action_bar,
            textvariable=self.confidence_var,
            command=self.show_confidence_details,
            width=36,
        )
        self.confidence_button.pack(side="right")
        self.confidence_button.state(["disabled"])
        ttk.Button(action_bar, text="Открыть CSV", command=self.open_file).pack(side="left", padx=(0, 10))
        self.prev_button = ttk.Button(action_bar, text="←", width=4, command=self.prev_batch)
        self.prev_button.pack(side="left", padx=(0, 4))
        self.next_button = ttk.Button(action_bar, text="→", width=4, command=self.next_batch)
        self.next_button.pack(side="left", padx=(0, 10))
        self.batch_results_button = ttk.Button(action_bar, text="Результаты Batch", command=self.open_batch_results_window)
        self.batch_results_button.pack(side="left", padx=(0, 10))
        ttk.Label(action_bar, textvariable=self.batch_var).pack(side="left", padx=(0, 14))

        info_bar = ttk.Frame(controls)
        info_bar.pack(fill="x", pady=(6, 0))
        ttk.Label(info_bar, textvariable=self.file_var).pack(side="left", padx=(0, 20))
        ttk.Label(info_bar, textvariable=self.omega_var).pack(side="left")
        ttk.Label(info_bar, textvariable=self.integration_var).pack(side="left", padx=(20, 0))
        ttk.Label(info_bar, textvariable=self.gamma_var).pack(side="left", padx=(20, 0))

        content = ttk.Frame(self.root, padding=(10, 0, 10, 0))
        content.pack(fill="both", expand=True)

        self.main_pane = ttk.Panedwindow(content, orient="horizontal")
        self.main_pane.pack(fill="both", expand=True)
        plot_frame = ttk.Frame(self.main_pane)
        sidebar = ttk.Frame(self.main_pane, padding=(10, 0, 0, 0))
        self.main_pane.add(plot_frame, weight=4)
        self.main_pane.add(sidebar, weight=1)

        # FigureCanvasTkAgg propagates the Figure's requested pixel size to Tk.
        # A fixed 1200x800 canvas made the whole UI wider than smaller screens.
        # These are only Tk's requested dimensions; the packed canvas expands
        # with the window.  Keep the request deliberately compact so panes can
        # also shrink cleanly on laptops and scaled Windows desktops.
        figure_width_px = max(650, int(self.initial_window_width * (0.52 if self.compact_ui else 0.62)))
        figure_height_px = max(460, int(self.initial_window_height - (240 if self.compact_ui else 180)))
        self.figure = Figure(figsize=(figure_width_px / 100.0, figure_height_px / 100.0), dpi=100)
        self.figure.subplots_adjust(left=0.055, right=0.985, top=0.955, bottom=0.065)
        grid = self.figure.add_gridspec(3, 2, height_ratios=[2.5, 1.0, 1.0], hspace=0.26, wspace=0.18)
        self.ax = self.figure.add_subplot(grid[0, :])
        self.preview_axes = []
        self.preview_specs = PREVIEW_WINDOWS[:]
        for index, spec in enumerate(self.preview_specs):
            row = 1 + index // 2
            col = index % 2
            preview_ax = self.figure.add_subplot(grid[row, col])
            self.preview_axes.append(preview_ax)

        self.canvas = FigureCanvasTkAgg(self.figure, master=plot_frame)
        self.canvas.mpl_connect("button_press_event", self.handle_manual_boundary_press)
        self.canvas.mpl_connect("motion_notify_event", self.handle_manual_boundary_motion)
        self.canvas.mpl_connect("button_release_event", self.handle_manual_boundary_release)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)
        toolbar_host = tk.Frame(plot_frame)
        toolbar_host.pack(fill="x")
        toolbar = NavigationToolbar2Tk(self.canvas, toolbar_host, pack_toolbar=False)
        toolbar.update()
        toolbar.pack(side="left", fill="x")

        self.sidebar_pane = ttk.Panedwindow(sidebar, orient="vertical")
        self.sidebar_pane.pack(fill="both", expand=True)
        batch_frame = ttk.LabelFrame(self.sidebar_pane, text="Batch", padding=(8, 8))
        self.sidebar_pane.add(batch_frame, weight=1)
        batch_columns = ("sample_name", "omega_value", "confidence")
        batch_tree_frame = ttk.Frame(batch_frame)
        batch_tree_frame.pack(fill="both", expand=True)
        self.batch_tree = ttk.Treeview(
            batch_tree_frame,
            columns=batch_columns,
            show="headings",
            height=8 if self.compact_ui else 14,
            selectmode="browse",
        )
        self.batch_tree.heading("sample_name", text="Образец")
        self.batch_tree.heading("omega_value", text="Omega-3")
        self.batch_tree.heading("confidence", text="Решение")
        batch_widths = (190, 70, 165) if self.compact_ui else (230, 80, 205)
        self.batch_tree.column("sample_name", width=batch_widths[0], minwidth=110, anchor="w", stretch=True)
        self.batch_tree.column("omega_value", width=batch_widths[1], minwidth=60, anchor="center", stretch=False)
        self.batch_tree.column("confidence", width=batch_widths[2], minwidth=125, anchor="w", stretch=True)
        self.batch_tree.pack(side="left", fill="both", expand=True)
        batch_scroll = ttk.Scrollbar(batch_tree_frame, orient="vertical", command=self.batch_tree.yview)
        batch_scroll.pack(side="right", fill="y")
        self.batch_tree.configure(yscrollcommand=batch_scroll.set)
        self.batch_tree.bind("<<TreeviewSelect>>", self.handle_batch_tree_selection)
        self._configure_quality_tags(self.batch_tree)

        details_frame = ttk.Frame(self.sidebar_pane)
        self.sidebar_pane.add(details_frame, weight=3)
        table_frame = ttk.LabelFrame(details_frame, text="Пики", padding=(8, 8))
        table_frame.pack(fill="both", expand=True, pady=(10, 0))
        # Keep the values used for manual review in the visible part of the
        # sidebar.  Area used to sit beyond the right edge of the table.
        cols = ["display_name", "area", "percent_area", "code", "expected_rt", "found_rt", "status"]
        peaks_tree_frame = ttk.Frame(table_frame)
        peaks_tree_frame.pack(fill="both", expand=True)
        self.tree = ttk.Treeview(
            peaks_tree_frame,
            columns=cols,
            show="headings",
            height=10 if self.compact_ui else 18,
        )
        headings = {
            "display_name": "Пик",
            "area": "Площадь",
            "percent_area": "%",
            "code": "Код",
            "expected_rt": "RT ожид.",
            "found_rt": "RT найден",
            "status": "Комментарий",
        }
        widths = {
            "display_name": 190,
            "area": 105,
            "percent_area": 58,
            "code": 85,
            "expected_rt": 82,
            "found_rt": 82,
            "status": 220,
        }
        if self.compact_ui:
            widths.update({
                "display_name": 150,
                "area": 85,
                "percent_area": 48,
                "code": 70,
                "expected_rt": 68,
                "found_rt": 68,
                "status": 150,
            })
        for c in cols:
            self.tree.heading(c, text=headings[c])
            width = widths[c]
            anchor = "w" if c in {"display_name", "status"} else "center"
            self.tree.column(c, width=width, anchor=anchor)
        self.tree.grid(row=0, column=0, sticky="nsew")
        peaks_scroll = ttk.Scrollbar(peaks_tree_frame, orient="vertical", command=self.tree.yview)
        peaks_scroll.grid(row=0, column=1, sticky="ns")
        peaks_scroll_x = ttk.Scrollbar(peaks_tree_frame, orient="horizontal", command=self.tree.xview)
        peaks_scroll_x.grid(row=1, column=0, sticky="ew")
        peaks_tree_frame.rowconfigure(0, weight=1)
        peaks_tree_frame.columnconfigure(0, weight=1)
        self.tree.configure(yscrollcommand=peaks_scroll.set, xscrollcommand=peaks_scroll_x.set)
        self.tree.bind("<<TreeviewSelect>>", self.handle_target_selection)

        manual_frame = ttk.LabelFrame(details_frame, text="Ручная интеграция", padding=(8, 8))
        manual_frame.pack(fill="x", pady=(10, 0))
        ttk.Label(manual_frame, text="Start RT").grid(row=0, column=0, sticky="w")
        ttk.Entry(manual_frame, textvariable=self.manual_start_var, width=12).grid(row=0, column=1, sticky="ew", padx=(6, 10))
        ttk.Label(manual_frame, text="End RT").grid(row=0, column=2, sticky="w")
        ttk.Entry(manual_frame, textvariable=self.manual_end_var, width=12).grid(row=0, column=3, sticky="ew", padx=(6, 0))
        ttk.Button(manual_frame, text="Взять текущие", command=self.load_selected_integration_bounds).grid(row=1, column=0, columnspan=2, sticky="ew", pady=(8, 0), padx=(0, 6))
        ttk.Button(manual_frame, text="Применить", command=self.apply_manual_integration).grid(row=1, column=2, columnspan=2, sticky="ew", pady=(8, 0))
        manual_frame.columnconfigure(1, weight=1)
        manual_frame.columnconfigure(3, weight=1)

        status = ttk.Label(self.root, textvariable=self.status_var, anchor="w", padding=(10, 4))
        status.pack(fill="x")
        self.update_batch_navigation()
        self.root.after_idle(self._apply_responsive_pane_positions)

    def _apply_responsive_pane_positions(self):
        """Give the plot and both sidebar sections useful space on this screen."""
        try:
            self.root.update_idletasks()
            pane_width = self.main_pane.winfo_width()
            pane_height = self.sidebar_pane.winfo_height()
            if pane_width > 1:
                plot_fraction = 0.68 if self.compact_ui else 0.72
                self.main_pane.sashpos(0, max(620, int(pane_width * plot_fraction)))
            if pane_height > 1:
                self.sidebar_pane.sashpos(0, max(150, int(pane_height * 0.30)))
        except tk.TclError:
            pass

    def update_batch_navigation(self):
        total = len(self.loaded_batches)
        if total <= 0:
            self.batch_var.set("Series: —")
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

        self.batch_var.set(f"Series: {self.current_batch_index + 1}/{total}")
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
        if risk_score >= 95:
            recommendation = "СТОП: откройте отмеченные пики и переинтегрируйте их вручную"
        elif risk_score >= 85:
            recommendation = "ПРОВЕРИТЬ: откройте отмеченные пики; при неверной границе переинтегрируйте"
        else:
            recommendation = "ГОТОВО: по проверкам судьи ручная интеграция не требуется"
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

    def handle_target_selection(self, event=None):
        selection = self.tree.selection()
        self.selected_target_code = selection[0] if selection else None
        self.load_selected_integration_bounds(silent=True)
        if self.df_processed is not None:
            self.update_plot()
            return

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
        self.manual_start_var.set("" if not np.isfinite(start_x) else f"{float(start_x):.5f}")
        self.manual_end_var.set("" if not np.isfinite(end_x) else f"{float(end_x):.5f}")
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
            messagebox.showerror("Ручная интеграция", "Start RT и End RT должны быть числами.", parent=self.root)
            return
        if not (np.isfinite(start_x) and np.isfinite(end_x) and end_x > start_x):
            messagebox.showerror("Ручная интеграция", "End RT должен быть больше Start RT.", parent=self.root)
            return

        x_col = _get_x_column_name(self.df_processed)
        x = self.df_processed[x_col].to_numpy(dtype=float)
        y = self.df_processed["y_corrected"].to_numpy(dtype=float)
        if start_x < float(np.nanmin(x)) or end_x > float(np.nanmax(x)):
            messagebox.showerror("Ручная интеграция", "Границы вне диапазона текущей хроматограммы.", parent=self.root)
            return

        start_idx = int(np.searchsorted(x, start_x, side="left"))
        end_idx = int(np.searchsorted(x, end_x, side="right") - 1)
        start_idx = max(0, min(start_idx, len(x) - 2))
        end_idx = max(start_idx + 1, min(end_idx, len(x) - 1))
        segment_y = np.clip(y[start_idx:end_idx + 1], 0.0, None)
        segment_x = x[start_idx:end_idx + 1]
        area = float(np.trapezoid(segment_y, segment_x))
        apex_idx = int(start_idx + np.argmax(segment_y))

        row_mask = self.matched_targets_df["code"] == self.selected_target_code
        if not row_mask.any():
            messagebox.showerror("Ручная интеграция", f"Пик {self.selected_target_code} не найден в таблице.", parent=self.root)
            return
        row_idx = self.matched_targets_df.index[row_mask][0]
        old_status = str(self.matched_targets_df.at[row_idx, "status"] or "")
        manual_status = old_status if "manual" in old_status else f"{old_status}_manual" if old_status else "manual"
        self.matched_targets_df.at[row_idx, "integration_start_x"] = float(x[start_idx])
        self.matched_targets_df.at[row_idx, "integration_end_x"] = float(x[end_idx])
        self.matched_targets_df.at[row_idx, "found_rt"] = float(x[apex_idx])
        self.matched_targets_df.at[row_idx, "area"] = area
        self.matched_targets_df.at[row_idx, "matched_peak_id"] = np.nan
        self.matched_targets_df.at[row_idx, "match_score"] = np.nan
        self.matched_targets_df.at[row_idx, "status"] = manual_status
        self.matched_targets_df = _recompute_matched_percent_area(self.matched_targets_df)
        self.manual_start_var.set(f"{float(x[start_idx]):.5f}")
        self.manual_end_var.set(f"{float(x[end_idx]):.5f}")
        self.selected_target_code = str(self.selected_target_code)
        # Keep the existing axes intact: clearing/recreating them also resets
        # Matplotlib toolbar zoom history even when xlim/ylim are restored.
        self.refresh_peaks(redraw_plot=False)
        if self._manual_overlay_artists:
            self._manual_drag_axis = None
            self._manual_drag_pending_bounds = (float(x[start_idx]), float(x[end_idx]))
            self._flush_manual_drag_overlay()
        else:
            self.update_plot(preserve_view=True)
        if self.selected_target_code in self.tree.get_children():
            self.tree.selection_set(self.selected_target_code)
            self.tree.focus(self.selected_target_code)
        self.status_var.set(
            f"Ручная интеграция {self.selected_target_code}: {x[start_idx]:.5f}–{x[end_idx]:.5f}, area {area:.4f}"
        )

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

    def handle_manual_boundary_press(self, event):
        if event.button != 1 or event.inaxes not in [self.ax, *getattr(self, "preview_axes", [])]:
            return
        if event.xdata is None or self.df_processed is None or not self.selected_target_code:
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
            self.manual_start_var.set(f"{x_value:.5f}")
        else:
            if np.isfinite(start_x):
                x_value = max(x_value, float(start_x) + 1e-5)
            self.manual_end_var.set(f"{x_value:.5f}")
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
                artists["fill"] = axis.fill_between(
                    x[start_idx:end_idx + 1],
                    0.0,
                    fill_y[start_idx:end_idx + 1],
                    color="#ff4d6d",
                    alpha=0.36 if artists["compact"] else 0.40,
                    linewidth=0.0,
                    zorder=3,
                )
                visual_fill_artist = artists.get("visual_fill")
                if visual_fill_artist is not None:
                    visual_fill_artist.remove()
                    selected_match = self.matched_targets_df[
                        self.matched_targets_df["code"] == self.selected_target_code
                    ]
                    if not selected_match.empty:
                        visual_row = selected_match.iloc[0].copy()
                        visual_row["integration_start_x"] = start_x
                        visual_row["integration_end_x"] = end_x
                        visual_start_x, visual_end_x = self._visual_peak_footprint_bounds(
                            visual_row,
                            x,
                            artists["y_smooth"],
                        )
                        visual_start_idx = int(np.argmin(np.abs(x - visual_start_x)))
                        visual_end_idx = int(np.argmin(np.abs(x - visual_end_x)))
                        artists["visual_fill"] = axis.fill_between(
                            x[visual_start_idx:visual_end_idx + 1],
                            0.0,
                            artists["visual_fill_y"][visual_start_idx:visual_end_idx + 1],
                            color="#ff7c96",
                            alpha=0.26 if artists["compact"] else 0.28,
                            linewidth=0.0,
                            zorder=2.7,
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
        if self.current_batch_index <= 0:
            return
        self.load_batch_at_index(self.current_batch_index - 1)

    def next_batch(self):
        if self.current_batch_index >= len(self.loaded_batches) - 1:
            return
        self.load_batch_at_index(self.current_batch_index + 1)

    def process_batch(self, batch: dict):
        if batch.get("processed_df") is not None:
            return batch
        batch.update(process_chromatogram_batch(batch["dataframe"], self.reference_targets))
        return batch

    def load_batch_at_index(self, index: int):
        if index < 0 or index >= len(self.loaded_batches):
            return

        batch = self.loaded_batches[index]
        self.process_batch(batch)
        self.current_batch_index = index
        self.current_sample_name = batch["sample_name"]
        self.file_var.set(
            f"Файл: {self.current_file.name} | Sample: {self.current_sample_name} | Source: {batch.get('file_name', '')}"
        )

        self.df_processed = batch["processed_df"]
        self.best_window = batch["best_window"]
        self.peaks_df = batch["peaks_df"]
        self.matched_targets_df = batch["matched_targets_df"]
        self.current_rt_shift = batch["rt_shift"]
        available_codes = set(self.matched_targets_df.get("code", pd.Series(dtype=str)).astype(str))
        if self.selected_target_code in available_codes:
            self.load_selected_integration_bounds(silent=True)
        else:
            self.selected_target_code = None
            self.manual_start_var.set("")
            self.manual_end_var.set("")
        self.update_batch_navigation()
        self.refresh_peaks()

    @staticmethod
    def _configure_quality_tags(tree: ttk.Treeview):
        tree.tag_configure("quality_stop", background="#ffd9dc", foreground="#8a1020")
        tree.tag_configure("quality_check", background="#fff0c7", foreground="#6b4a00")
        tree.tag_configure("quality_good", background="#def3e5", foreground="#155b32")
        tree.tag_configure("quality_pending", foreground="#6f7780")

    @staticmethod
    def _operator_quality(confidence: dict | None) -> tuple[str, str]:
        """Translate the judge into one unambiguous operator action."""
        if not isinstance(confidence, dict) or not confidence:
            return "ОЖИДАНИЕ — ещё не рассчитано", "quality_pending"
        risk = confidence.get("high_error_risk", {})
        risk_score = risk.get("score", 0) if isinstance(risk, dict) else 0
        retry = confidence.get("structural_retry", {})
        if risk_score >= 95:
            return "СТОП — переинтегрировать", "quality_stop"
        if risk_score >= 85:
            risky_codes = ", ".join(risk.get("peak_codes", [])) if isinstance(risk, dict) else ""
            suffix = f": {risky_codes}" if risky_codes else " отмеченные пики"
            return f"ПРОВЕРИТЬ{suffix}", "quality_check"
        if isinstance(retry, dict) and retry.get("accepted"):
            return "ГОТОВО — перепроверено автоматически", "quality_good"
        return "ГОТОВО — правка не нужна", "quality_good"

    def build_batch_results_rows(self, process_all: bool = False):
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
        for index, sample_name, value_text, confidence_text, quality_tag in self.build_batch_results_rows(process_all=process_all):
            display_name = self._sample_number(sample_name) if tree is self.batch_results_tree else sample_name
            values = (display_name, value_text, confidence_text) if show_confidence else (display_name, value_text)
            tree.insert("", "end", iid=str(index), values=values, tags=(quality_tag,))
        if selected_iid is not None and tree.exists(selected_iid):
            self._batch_tree_syncing = True
            tree.selection_set(selected_iid)
            tree.focus(selected_iid)
            tree.see(selected_iid)
            self._batch_tree_syncing = False

    def populate_main_batch_tree(self):
        self._populate_batch_tree_widget(self.batch_tree, process_all=False)

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
        self.batch_results_window.title("Batch Results")
        self.batch_results_window.geometry("860x600")
        self.batch_results_window.minsize(650, 420)

        frame = ttk.Frame(self.batch_results_window, padding=10)
        frame.pack(fill="both", expand=True)

        toolbar = ttk.Frame(frame)
        toolbar.pack(fill="x", pady=(0, 8))
        ttk.Button(toolbar, text="Copy Selected", command=lambda: self.copy_batch_results(selected_only=True)).pack(side="left")
        ttk.Button(toolbar, text="Copy All", command=lambda: self.copy_batch_results(selected_only=False)).pack(side="left", padx=(8, 0))
        ttk.Button(toolbar, text="Open Selected", command=self.jump_to_batch_from_results).pack(side="left", padx=(8, 0))

        tree_frame = ttk.Frame(frame)
        tree_frame.pack(fill="both", expand=True)
        columns = ("sample_name", "omega_value", "confidence")
        self.batch_results_tree = ttk.Treeview(tree_frame, columns=columns, show="headings", height=18, selectmode="extended")
        self.batch_results_tree.heading("sample_name", text="Номер образца")
        self.batch_results_tree.heading("omega_value", text="Значение")
        self.batch_results_tree.heading("confidence", text="Проверка")
        self.batch_results_tree.column("sample_name", width=260, minwidth=140, anchor="w", stretch=True)
        self.batch_results_tree.column("omega_value", width=110, minwidth=80, anchor="center", stretch=False)
        self.batch_results_tree.column("confidence", width=430, minwidth=260, anchor="w", stretch=True)
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

    def open_file(self):
        file_path = filedialog.askopenfilename(title="Выберите CSV", filetypes=[("CSV", "*.csv *.CSV"), ("All", "*.*")])
        if not file_path:
            return
        try:
            self.current_file = Path(file_path)
            self.loaded_batches = omega_core.load_batches(self.current_file, cutoff_minutes=4.0)
            if self._preload_after_id is not None:
                self.root.after_cancel(self._preload_after_id)
                self._preload_after_id = None
            self._preload_batch_index = 0
            self.show_batch_progress_window()
            self._preload_after_id = self.root.after(50, self.preload_loaded_batches)
            self.status_var.set(
                f"Загружено проб: {len(self.loaded_batches)}. "
                "Запускаю последовательный расчёт."
            )
        except Exception as e:
            self.close_batch_progress_window()
            messagebox.showerror("Ошибка", str(e), parent=self.root)

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
            return

        index = self._preload_batch_index
        batch = self.loaded_batches[index]
        sample_name = batch.get("sample_name", batch.get("file_name", f"Проба {index + 1}"))
        self.batch_progress_label_var.set(f"Проба {index + 1} из {total}")
        self.batch_progress_detail_var.set(f"Сейчас анализируется: {sample_name}")
        if self.batch_progress_bar is not None:
            self.batch_progress_bar["maximum"] = total
            self.batch_progress_bar["value"] = index
        self.root.update_idletasks()
        self.status_var.set(f"Расчёт пробы: {index + 1}/{total}")
        try:
            self.process_batch(batch)
        except Exception as exc:
            self.close_batch_progress_window()
            self.status_var.set(f"Ошибка при расчёте пробы {index + 1}/{total}")
            messagebox.showerror("Ошибка анализа", f"{sample_name}\n\n{exc}", parent=self.root)
            return
        self._preload_batch_index += 1
        self._preload_after_id = self.root.after(25, self.preload_loaded_batches)

    def refresh_peaks(self, preserve_plot_view: bool = False, redraw_plot: bool = True):
        if self.df_processed is None:
            return
        current_batch = self.loaded_batches[self.current_batch_index] if self.loaded_batches else None
        omega = core_metrics.compute_omega(self.matched_targets_df)
        baseline_mode = current_batch.get("baseline_mode", "chebyshev") if current_batch is not None else "chebyshev"
        cluster_quality_score = core_metrics.compute_cluster_quality(self.matched_targets_df)
        confidence = core_metrics.assess_confidence(
            self.matched_targets_df,
            self.peaks_df,
            omega,
            baseline_mode,
            cluster_quality_score,
        )
        report_value = omega["omega3_trio"]
        if current_batch is not None:
            current_batch["processed_df"] = self.df_processed
            current_batch["best_window"] = self.best_window
            current_batch["peaks_df"] = self.peaks_df
            current_batch["matched_targets_df"] = self.matched_targets_df
            current_batch["rt_shift"] = self.current_rt_shift
            current_batch["omega"] = omega
            current_batch["omega_report"] = omega["omega3_trio"]
            current_batch["cluster_quality_score"] = cluster_quality_score
            current_batch["confidence"] = confidence
            report_value = current_batch["omega_report"]

        if np.isfinite(report_value):
            self.omega_var.set(
                f"Omega-3: {report_value:.2f}% | strict: {omega['omega3_trio_strict']:.2f}%"
            )
        else:
            self.omega_var.set("Omega-3: —")

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

        if redraw_plot:
            self.update_plot(preserve_view=preserve_plot_view)
        self.update_table()

        self.status_var.set(
            f"RT shift: {self.current_rt_shift:+.3f} min | matched {int(self.matched_targets_df['matched_peak_id'].notna().sum())}/{len(self.reference_targets)}"
        )
        self.integration_var.set(
            f"Integration: {len(self.peaks_df)} peaks | SG {self.best_window}"
        )
        self.populate_main_batch_tree()
        if self.batch_results_window is not None and self.batch_results_window.winfo_exists():
            self.populate_batch_results_tree()

    def _resolve_selected_plot_items(self):
        selected_row = None
        selected_peak = None
        if self.selected_target_code and not self.matched_targets_df.empty:
            selected_match = self.matched_targets_df[self.matched_targets_df["code"] == self.selected_target_code]
            if not selected_match.empty:
                selected_row = selected_match.iloc[0]
                matched_peak_id = pd.to_numeric(selected_match["matched_peak_id"], errors="coerce").iloc[0]
                if np.isfinite(matched_peak_id) and not self.peaks_df.empty:
                    selected_peak_match = self.peaks_df[self.peaks_df["peak_id"] == int(matched_peak_id)]
                    if not selected_peak_match.empty:
                        selected_peak = selected_peak_match.iloc[0]
        return selected_row, selected_peak

    def _visual_peak_footprint_bounds(
        self,
        target_row,
        x: np.ndarray,
        y_smooth: np.ndarray,
    ) -> tuple[float, float]:
        """Extend only the painted footprint to the visible feet of a peak.

        The calculation continues to use integration_start/end_x.  This helper
        prevents a non-zero first sample from looking like a peak cut by a
        vertical knife while keeping the validated numeric boundaries intact.
        """
        start_x = pd.to_numeric(pd.Series([target_row.get("integration_start_x")]), errors="coerce").iloc[0]
        end_x = pd.to_numeric(pd.Series([target_row.get("integration_end_x")]), errors="coerce").iloc[0]
        apex_x = pd.to_numeric(pd.Series([target_row.get("found_rt")]), errors="coerce").iloc[0]
        if not np.all(np.isfinite([start_x, end_x, apex_x])) or end_x <= start_x:
            return float(start_x), float(end_x)

        start_idx = int(np.argmin(np.abs(x - float(start_x))))
        end_idx = int(np.argmin(np.abs(x - float(end_x))))
        apex_idx = int(np.argmin(np.abs(x - float(apex_x))))
        if not (start_idx < apex_idx < end_idx):
            return float(start_x), float(end_x)

        ordered_rts = np.sort(pd.to_numeric(
            self.matched_targets_df.get("found_rt", pd.Series(dtype=float)),
            errors="coerce",
        ).dropna().to_numpy(dtype=float))
        previous = ordered_rts[ordered_rts < float(apex_x) - 1e-9]
        following = ordered_rts[ordered_rts > float(apex_x) + 1e-9]
        left_limit_x = float(apex_x) - 0.040
        right_limit_x = float(apex_x) + 0.040
        if previous.size and float(apex_x) - float(previous[-1]) <= 0.12:
            left_limit_x = max(left_limit_x, 0.5 * (float(previous[-1]) + float(apex_x)))
        if following.size and float(following[0]) - float(apex_x) <= 0.12:
            right_limit_x = min(right_limit_x, 0.5 * (float(following[0]) + float(apex_x)))
        left_limit_x = max(left_limit_x, float(start_x) - 0.024)
        right_limit_x = min(right_limit_x, float(end_x) + 0.024)
        left_limit_idx = max(0, int(np.searchsorted(x, left_limit_x, side="left")))
        right_limit_idx = min(len(x) - 1, int(np.searchsorted(x, right_limit_x, side="right") - 1))

        local_left = max(0, left_limit_idx - 8)
        local_right = min(len(x), right_limit_idx + 9)
        local_signal = np.asarray(y_smooth[local_left:local_right], dtype=float)
        finite_local = local_signal[np.isfinite(local_signal)]
        noise = float(np.median(np.abs(finite_local - np.median(finite_local))) * 1.4826) if finite_local.size else 0.0
        apex_height = max(float(y_smooth[apex_idx]), 1.0)
        visible_step = max(noise * 0.35, abs(apex_height) * 0.004, 1e-9)

        visual_start_idx = start_idx
        if left_limit_idx < start_idx:
            segment = np.asarray(y_smooth[left_limit_idx:start_idx + 1], dtype=float)
            if np.isfinite(segment).any():
                candidate = left_limit_idx + int(np.nanargmin(segment))
                if candidate < start_idx and float(y_smooth[start_idx] - y_smooth[candidate]) > visible_step:
                    visual_start_idx = candidate

        visual_end_idx = end_idx
        if end_idx < right_limit_idx:
            segment = np.asarray(y_smooth[end_idx:right_limit_idx + 1], dtype=float)
            if np.isfinite(segment).any():
                candidate = end_idx + int(np.nanargmin(segment))
                if candidate > end_idx and float(y_smooth[end_idx] - y_smooth[candidate]) > visible_step:
                    visual_end_idx = candidate

        return float(x[visual_start_idx]), float(x[visual_end_idx])

    def _draw_chromatogram_axis(
        self,
        axis,
        x: np.ndarray,
        y: np.ndarray,
        y_smooth: np.ndarray,
        fill_y: np.ndarray,
        marker_y: np.ndarray,
        selected_row,
        selected_peak,
        title: str,
        x_min=None,
        x_max=None,
        compact: bool = False,
        normalized: bool = False,
    ):
        axis.clear()
        axis.set_facecolor("#fcfcfc")
        y_draw = y
        y_smooth_draw = y_smooth
        fill_y_draw = fill_y
        marker_y_draw = marker_y
        if normalized:
            visible_mask = np.ones(len(x), dtype=bool)
            if x_min is not None:
                visible_mask &= x >= float(x_min)
            if x_max is not None:
                visible_mask &= x <= float(x_max)
            if np.any(visible_mask):
                local_candidates = [
                    np.abs(y_smooth[visible_mask]),
                    np.abs(marker_y[visible_mask]),
                    np.abs(fill_y[visible_mask]),
                ]
                local_scale = max(
                    1e-9,
                    max(float(np.nanmax(values)) for values in local_candidates if values.size > 0),
                )
                y_draw = y / local_scale
                y_smooth_draw = y_smooth / local_scale
                fill_y_draw = fill_y / local_scale
                marker_y_draw = marker_y / local_scale
        axis.axhline(0.0, color="#777777", linewidth=0.8, alpha=0.55)
        axis.grid(color="#d9d9d9", linewidth=0.45, alpha=0.55)
        axis.plot(x, y_draw, linewidth=0.95 if compact else 1.0, color="#2a5b84", alpha=0.50, label="Corrected")
        axis.plot(x, y_smooth_draw, linewidth=1.05 if compact else 1.2, color="#111111", alpha=0.92, label="Smoothed")

        if not self.peaks_df.empty:
            for _, peak in self.peaks_df.iterrows():
                peak_start_x = float(peak["start_x"])
                peak_end_x = float(peak["end_x"])
                peak_apex_x = float(peak["apex_x"])
                if x_min is not None and peak_end_x < x_min:
                    continue
                if x_max is not None and peak_start_x > x_max:
                    continue
                start_idx = int(peak["start_idx"])
                end_idx = int(peak["end_idx"])
                apex_idx = int(peak["apex_idx"])
                axis.axvline(peak_start_x, color="#caa25a", linewidth=0.45 if compact else 0.5, alpha=0.18)
                axis.axvline(peak_end_x, color="#caa25a", linewidth=0.45 if compact else 0.5, alpha=0.18)
                axis.scatter(x[apex_idx], marker_y_draw[apex_idx], s=10 if compact else 16, color="#b84a35", alpha=0.75, zorder=4)

        if not self.matched_targets_df.empty:
            for _, target_row in self.matched_targets_df.iterrows():
                start_x = pd.to_numeric(pd.Series([target_row.get("integration_start_x")]), errors="coerce").iloc[0]
                end_x = pd.to_numeric(pd.Series([target_row.get("integration_end_x")]), errors="coerce").iloc[0]
                if not (np.isfinite(start_x) and np.isfinite(end_x)):
                    continue
                if x_min is not None and end_x < x_min:
                    continue
                if x_max is not None and start_x > x_max:
                    continue
                start_idx = int(np.argmin(np.abs(x - float(start_x))))
                end_idx = int(np.argmin(np.abs(x - float(end_x))))
                if end_idx <= start_idx:
                    continue
                visual_start_x, visual_end_x = self._visual_peak_footprint_bounds(target_row, x, y_smooth)
                visual_start_idx = int(np.argmin(np.abs(x - visual_start_x)))
                visual_end_idx = int(np.argmin(np.abs(x - visual_end_x)))
                visual_fill = np.clip(y_smooth_draw, 0.0, None)
                if visual_end_idx > visual_start_idx:
                    axis.fill_between(
                        x[visual_start_idx:visual_end_idx + 1],
                        0.0,
                        visual_fill[visual_start_idx:visual_end_idx + 1],
                        color="#f2b134",
                        alpha=0.18 if compact else 0.20,
                        linewidth=0.0,
                        zorder=1.8,
                    )
                axis.fill_between(
                    x[start_idx:end_idx + 1],
                    0.0,
                    fill_y_draw[start_idx:end_idx + 1],
                    color="#f2b134",
                    alpha=0.22 if compact else 0.24,
                    linewidth=0.0,
                    zorder=2,
                )
                axis.axvline(float(start_x), color="#d6a033", linewidth=0.65 if compact else 0.75, alpha=0.32, zorder=3)
                axis.axvline(float(end_x), color="#d6a033", linewidth=0.65 if compact else 0.75, alpha=0.32, zorder=3)

        if selected_row is not None:
            start_x = pd.to_numeric(pd.Series([selected_row.get("integration_start_x") if selected_row is not None else np.nan]), errors="coerce").iloc[0]
            end_x = pd.to_numeric(pd.Series([selected_row.get("integration_end_x") if selected_row is not None else np.nan]), errors="coerce").iloc[0]
            manual_start_x, manual_end_x = self._manual_bounds_from_vars()
            if np.isfinite(manual_start_x) and np.isfinite(manual_end_x):
                start_x, end_x = manual_start_x, manual_end_x
            apex_x = pd.to_numeric(pd.Series([selected_row.get("found_rt") if selected_row is not None else np.nan]), errors="coerce").iloc[0]
            if (not np.isfinite(start_x) or not np.isfinite(end_x)) and selected_peak is not None:
                start_x = float(selected_peak["start_x"])
                end_x = float(selected_peak["end_x"])
            if not np.isfinite(apex_x) and selected_peak is not None:
                apex_x = float(selected_peak["apex_x"])
            if not np.isfinite(apex_x) and np.isfinite(start_x) and np.isfinite(end_x):
                apex_x = 0.5 * (float(start_x) + float(end_x))
            if np.isfinite(start_x) and np.isfinite(end_x) and np.isfinite(apex_x) and (x_min is None or end_x >= x_min) and (x_max is None or start_x <= x_max):
                start_idx = int(np.argmin(np.abs(x - float(start_x))))
                end_idx = int(np.argmin(np.abs(x - float(end_x))))
                apex_idx = int(np.argmin(np.abs(x - float(apex_x))))
                visual_start_x, visual_end_x = self._visual_peak_footprint_bounds(selected_row, x, y_smooth)
                visual_start_idx = int(np.argmin(np.abs(x - visual_start_x)))
                visual_end_idx = int(np.argmin(np.abs(x - visual_end_x)))
                visual_fill = np.clip(y_smooth_draw, 0.0, None)
                selected_visual_fill_artist = None
                if visual_end_idx > visual_start_idx:
                    selected_visual_fill_artist = axis.fill_between(
                        x[visual_start_idx:visual_end_idx + 1],
                        0.0,
                        visual_fill[visual_start_idx:visual_end_idx + 1],
                        color="#ff7c96",
                        alpha=0.26 if compact else 0.28,
                        linewidth=0.0,
                        zorder=2.7,
                    )
                selected_fill_artist = axis.fill_between(
                    x[start_idx:end_idx + 1],
                    0.0,
                    fill_y_draw[start_idx:end_idx + 1],
                    color="#ff4d6d",
                    alpha=0.36 if compact else 0.40,
                    linewidth=0.0,
                    zorder=3,
                )
                start_line = axis.axvline(float(start_x), color="#ff4d6d", linewidth=1.0 if compact else 1.2, alpha=0.85, zorder=5)
                end_line = axis.axvline(float(end_x), color="#ff4d6d", linewidth=1.0 if compact else 1.2, alpha=0.85, zorder=5)
                self._manual_overlay_artists[axis] = {
                    "fill": selected_fill_artist,
                    "start_line": start_line,
                    "end_line": end_line,
                    "x": x,
                    "fill_y": fill_y_draw,
                    "y_smooth": y_smooth,
                    "visual_fill_y": visual_fill,
                    "visual_fill": selected_visual_fill_artist,
                    "marker_y": marker_y_draw,
                    "compact": compact,
                }
                selected_marker_artist = axis.scatter(
                    x[apex_idx],
                    marker_y_draw[apex_idx],
                    s=44 if compact else 70,
                    facecolor="#fff3f6",
                    edgecolor="#ff4d6d",
                    linewidth=1.3 if compact else 1.6,
                    zorder=6,
                )
                self._manual_overlay_artists[axis]["marker"] = selected_marker_artist

        visible_codes = []
        if not self.matched_targets_df.empty:
            labeled = self.matched_targets_df[self.matched_targets_df["matched_peak_id"].notna()].copy()
            for _, row in labeled.iterrows():
                found_rt = float(row["found_rt"])
                if x_min is not None and found_rt < x_min:
                    continue
                if x_max is not None and found_rt > x_max:
                    continue
                apex_idx = int(np.argmin(np.abs(x - found_rt)))
                visible_codes.append(str(row["code"]))
                axis.text(
                    found_rt,
                    marker_y_draw[apex_idx] + max(
                        np.nanmax(marker_y_draw) * 0.010,
                        0.06 if normalized else (50.0 if compact else 80.0),
                    ),
                    str(row["code"]),
                    fontsize=6.4 if compact else 7,
                    rotation=90,
                    ha="center",
                    va="bottom",
                    color="#2c3e50",
                    alpha=0.92,
                )

        if selected_row is not None and pd.notna(selected_row.get("found_rt")):
            selected_rt = float(selected_row["found_rt"])
            if (x_min is None or selected_rt >= x_min) and (x_max is None or selected_rt <= x_max):
                apex_idx = int(np.argmin(np.abs(x - selected_rt)))
                label_text = f"{selected_row['code']}  RT {selected_rt:.4f}"
                selected_annotation = axis.annotate(
                    label_text,
                    xy=(selected_rt, marker_y_draw[apex_idx]),
                    xytext=(10, 12 if compact else 14),
                    textcoords="offset points",
                    fontsize=7 if compact else 8,
                    color="#7a1028",
                    bbox={"boxstyle": "round,pad=0.25", "facecolor": "#fff3f6", "edgecolor": "#ff4d6d", "alpha": 0.95},
                    arrowprops={"arrowstyle": "->", "color": "#ff4d6d", "lw": 0.9},
                    zorder=7,
                )
                if axis in self._manual_overlay_artists:
                    self._manual_overlay_artists[axis]["annotation"] = selected_annotation

        if x_min is not None and x_max is not None:
            axis.set_xlim(float(x_min), float(x_max))
            if visible_codes:
                peak_text = ", ".join(visible_codes)
                axis.text(
                    0.01,
                    0.98,
                    peak_text,
                    transform=axis.transAxes,
                    ha="left",
                    va="top",
                    fontsize=7,
                    color="#304860",
                    bbox={"boxstyle": "round,pad=0.20", "facecolor": "#ffffff", "edgecolor": "#d5dde5", "alpha": 0.85},
                )
            if normalized:
                visible_mask = (x >= float(x_min)) & (x <= float(x_max))
                local_min = float(np.nanmin(y_draw[visible_mask])) if np.any(visible_mask) else 0.0
                local_max = float(np.nanmax(np.maximum.reduce([
                    np.asarray(y_draw[visible_mask]),
                    np.asarray(y_smooth_draw[visible_mask]),
                    np.asarray(fill_y_draw[visible_mask]),
                ]))) if np.any(visible_mask) else 1.0
                axis.set_ylim(min(-0.08, local_min * 1.08), max(1.05, local_max * 1.12))
        axis.set_title(title, fontsize=9 if compact else 11, pad=6)
        axis.tick_params(labelsize=7 if compact else 9)
        axis.set_xlabel("Time, min", fontsize=8 if compact else 10)
        axis.set_ylabel("Norm." if normalized else "Signal", fontsize=8 if compact else 10)
        for spine in axis.spines.values():
            spine.set_color("#b8c2cc")
            spine.set_linewidth(0.8)

    def update_plot(self, preserve_view: bool = False):
        axes = [self.ax, *getattr(self, "preview_axes", [])]
        saved_views = []
        if preserve_view:
            for axis in axes:
                saved_views.append((axis.get_xlim(), axis.get_ylim()))
        x_col = _get_x_column_name(self.df_processed)
        x = self.df_processed[x_col].to_numpy(dtype=float)
        y = self.df_processed["y_corrected"].to_numpy(dtype=float)
        y_smooth = self.df_processed["y_smooth"].to_numpy(dtype=float)
        fill_y = np.clip(y, 0.0, None)
        marker_y = y
        selected_row, selected_peak = self._resolve_selected_plot_items()
        self._manual_overlay_artists = {}

        self._draw_chromatogram_axis(
            self.ax,
            x=x,
            y=y,
            y_smooth=y_smooth,
            fill_y=fill_y,
            marker_y=marker_y,
            selected_row=selected_row,
            selected_peak=selected_peak,
            title="Общая хроматограмма",
            compact=False,
        )
        self.ax.legend(loc="upper right", fontsize=8)

        for preview_ax, (label, x_min, x_max) in zip(self.preview_axes, self.preview_specs):
            self._draw_chromatogram_axis(
                preview_ax,
                x=x,
                y=y,
                y_smooth=y_smooth,
                fill_y=fill_y,
                marker_y=marker_y,
                selected_row=selected_row,
                selected_peak=selected_peak,
                title=f"Фрагмент основного графика: {label}",
                x_min=x_min,
                x_max=x_max,
                compact=True,
                normalized=True,
            )
        if preserve_view and len(saved_views) == len(axes):
            for axis, (x_limits, y_limits) in zip(axes, saved_views):
                axis.set_xlim(*x_limits)
                axis.set_ylim(*y_limits)
        self.canvas.draw_idle()

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
                "" if pd.isna(row.get("expected_rt")) else f"{row['expected_rt']:.4f}",
                "" if pd.isna(row.get("found_rt")) else f"{row['found_rt']:.4f}",
                row.get("status", ""),
            )
            self.tree.insert("", "end", iid=code, values=vals)
        if self.selected_target_code not in available_codes:
            self.selected_target_code = None
        if self.selected_target_code is not None:
            self.tree.selection_set(self.selected_target_code)
            self.tree.focus(self.selected_target_code)


if __name__ == "__main__":
    root = tk.Tk()
    app = ChromatogramApp(root)
    root.mainloop()
