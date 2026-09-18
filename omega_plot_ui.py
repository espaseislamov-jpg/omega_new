"""Chromatogram rendering and viewport state; no analytical calculations."""
import numpy as np
import pandas as pd
from matplotlib.ticker import MaxNLocator, ScalarFormatter, AutoMinorLocator
from omega_core import integration


def _get_x_column_name(df: pd.DataFrame) -> str:
    if "x_corrected" in df.columns:
        return "x_corrected"
    if "x" in df.columns:
        return "x"
    raise KeyError("Expected an x or x_corrected column")


class ChromatogramPlotUI:
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
        self, target_row, x: np.ndarray, y_smooth: np.ndarray,
    ) -> tuple[float, float]:
        """Paint only the interval actually used for numerical integration."""
        start = pd.to_numeric(target_row.get("integration_start_x"), errors="coerce")
        end = pd.to_numeric(target_row.get("integration_end_x"), errors="coerce")
        return float(start), float(end)

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
        axis.set_axis_on()
        axis.set_facecolor("#FFFFFF")
        overview = axis is getattr(self, "ax", None)
        y_draw = y
        fill_y_draw = y
        marker_y_draw = marker_y
        if normalized:
            visible_mask = np.ones(len(x), dtype=bool)
            if x_min is not None:
                visible_mask &= x >= float(x_min)
            if x_max is not None:
                visible_mask &= x <= float(x_max)
            if np.any(visible_mask):
                local_candidates = [
                    np.abs(y[visible_mask]),
                    np.abs(marker_y[visible_mask]),
                    np.abs(fill_y[visible_mask]),
                ]
                local_scale = max(
                    1e-9,
                    max(float(np.nanmax(values)) for values in local_candidates if values.size > 0),
                )
                y_draw = y / local_scale
                fill_y_draw = y / local_scale
                marker_y_draw = marker_y / local_scale
        axis.axhline(0.0, color="#777777", linewidth=0.8, alpha=0.55)
        axis.grid(color="#E1E7EE", linewidth=0.5, alpha=0.65)
        axis.plot(x, y_draw, linewidth=0.95 if compact else 1.0, color="#111111", alpha=1.0, label="Corrected")

        if not self.peaks_df.empty:
            assigned_ids = set(self.matched_targets_df.get("matched_peak_id", pd.Series(dtype=float)).dropna())
            for _, peak in self.peaks_df.iterrows():
                if peak["peak_id"] in assigned_ids:
                    continue
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
                axis.scatter(x[apex_idx], marker_y_draw[apex_idx], s=7 if compact else 10, color="#777777", alpha=0.35, zorder=2)

        if not self.matched_targets_df.empty:
            for _, target_row in self.matched_targets_df.iterrows():
                color = {"C18:2N6C": "#d99e22", "C18:1N9C": "#268bd2", "C18:3N3": "#36965e",
                         "C20:5": "#268bd2", "C22:6": "#d99e22", "C22:5": "#268bd2", "C22:4": "#36965e"}.get(target_row["code"], "#d99e22")
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
                if selected_row is not None:
                    color = "#AEBFD1"
                fill_x, exact_y = integration.interval_points(x, y_draw, float(start_x), float(end_x))
                axis.fill_between(fill_x, 0.0, np.maximum(exact_y, 0.0), color=color,
                                  alpha=0.30, linewidth=0.0, zorder=2)
                if not overview:
                    axis.axvline(float(start_x), color=color, linewidth=0.65, alpha=0.45, zorder=3)
                    axis.axvline(float(end_x), color=color, linewidth=0.65, alpha=0.45, zorder=3)
                apex = target_row.get("found_rt", np.nan)
                if np.isfinite(apex):
                    axis.scatter([apex], [np.interp(apex, x, marker_y_draw)], color=color, s=12 if compact else 18, zorder=4)

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
                fill_x, exact_y = integration.interval_points(x, fill_y_draw, float(start_x), float(end_x))
                selected_fill_artist = axis.fill_between(
                    fill_x,
                    0.0,
                    np.maximum(exact_y, 0.0),
                    color="#2463C4",
                    alpha=0.36 if compact else 0.40,
                    linewidth=0.0,
                    zorder=3,
                )
                start_line = axis.axvline(float(start_x), ymax=.97, marker='s', markevery=[1], markersize=5,
                                         color="#2463C4", linewidth=1.0 if compact else 1.2, alpha=0.85, zorder=5)
                end_line = axis.axvline(float(end_x), ymax=.97, marker='s', markevery=[1], markersize=5,
                                       color="#2463C4", linewidth=1.0 if compact else 1.2, alpha=0.85, zorder=5)
                self._manual_overlay_artists[axis] = {
                    "fill": selected_fill_artist,
                    "start_line": start_line,
                    "end_line": end_line,
                    "x": x,
                    "fill_y": fill_y_draw,
                    "marker_y": marker_y_draw,
                    "compact": compact,
                }
                selected_marker_artist = axis.scatter(
                    x[apex_idx],
                    marker_y_draw[apex_idx],
                    s=44 if compact else 70,
                    facecolor="#EDF4FF",
                    edgecolor="#2463C4",
                    linewidth=1.3 if compact else 1.6,
                    zorder=6,
                )
                self._manual_overlay_artists[axis]["marker"] = selected_marker_artist

        visible_codes = []
        if not self.matched_targets_df.empty:
            labeled = self.matched_targets_df[self.matched_targets_df["matched_peak_id"].notna()].copy()
            for label_index, (_, row) in enumerate(labeled.iterrows()):
                found_rt = float(row["found_rt"])
                if x_min is not None and found_rt < x_min:
                    continue
                if x_max is not None and found_rt > x_max:
                    continue
                apex_idx = int(np.argmin(np.abs(x - found_rt)))
                visible_codes.append(str(row["code"]))
                if overview or (selected_row is not None and row['code'] == selected_row['code']):
                    continue
                if compact:
                    continue
                axis.annotate(str(row['code']), xy=(found_rt, marker_y_draw[apex_idx]),
                              xytext=(0, 12+14*(label_index % 2)), textcoords='offset points',
                              fontsize=8, ha='center', color='#596B80', annotation_clip=True)


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
                    color="#174B91",
                    bbox={"boxstyle": "round,pad=0.25", "facecolor": "#EDF4FF", "edgecolor": "#2463C4", "alpha": 0.95},
                    arrowprops={"arrowstyle": "->", "color": "#2463C4", "lw": 0.9},
                    zorder=7,
                )
                if axis in self._manual_overlay_artists:
                    self._manual_overlay_artists[axis]["annotation"] = selected_annotation

        if x_min is not None and x_max is not None:
            axis.set_xlim(float(x_min), float(x_max))
            if visible_codes and compact:
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
                    np.asarray(fill_y_draw[visible_mask]),
                ]))) if np.any(visible_mask) else 1.0
                axis.set_ylim(min(-0.08, local_min * 1.08), max(1.15, local_max * 1.32))
        axis.set_title(title, fontsize=9 if compact or overview else 11, pad=10, loc="left", color="#223047")
        axis.tick_params(labelsize=8, colors="#596B80")
        if overview:
            axis.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
            axis.yaxis.get_offset_text().set_fontsize(8)
        axis.set_xlabel("Время, мин", fontsize=8 if compact else 9)
        axis.xaxis.set_major_locator(MaxNLocator(nbins=5 if compact else 8, min_n_ticks=3))
        formatter = ScalarFormatter(useOffset=False)
        formatter.set_scientific(False)
        axis.xaxis.set_major_formatter(formatter)
        axis.xaxis.set_minor_locator(AutoMinorLocator(2))
        axis.tick_params(axis='x', labelbottom=True, bottom=True, which='both')
        axis.set_ylabel("Отн. сигнал" if normalized else "Сигнал", fontsize=8 if compact else 9)
        for spine in axis.spines.values():
            spine.set_color("#b8c2cc")
            spine.set_linewidth(0.8)

    @staticmethod
    def _capture_axes_views(axes):
        return [
            (tuple(axis.get_xlim()), tuple(axis.get_ylim()))
            for axis in axes
        ]

    @staticmethod
    def _restore_axes_views(axes, saved_views):
        if len(saved_views) != len(axes):
            return
        for axis, (x_limits, y_limits) in zip(axes, saved_views):
            axis.set_xlim(*x_limits)
            axis.set_ylim(*y_limits)

    def update_plot(self, preserve_view: bool = False):
        axes = [self.ax, *getattr(self, "preview_axes", [])]
        saved_views = self._capture_axes_views(axes) if preserve_view else []
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
            title="",
            compact=False,
        )
        self.ax.set_ylabel("")
        self.ax.set_xlabel("Время, мин", fontsize=8)
        if not self.full_overview_var.get():
            # Display only the region containing assigned fatty acids. This is
            # view scaling, not a change to the signal, identity or integration.
            valid = self.matched_targets_df.dropna(subset=['found_rt'])
            if not valid.empty:
                lo = float(valid.found_rt.min())-.3
                hi = float(valid.found_rt.max())+.3
                if 'integration_start_x' in valid:
                    lo = min(lo, float(valid.integration_start_x.min())-.05)
                    hi = max(hi, float(valid.integration_end_x.max())+.05)
                mask = (x >= lo) & (x <= hi)
                if mask.any():
                    self.ax.set_xlim(max(x[0],lo), min(x[-1],hi))
                    low, high = float(np.min(y[mask])), float(np.max(y[mask]))
                    span = max(high-low, 1.)
                    self.ax.set_ylim(min(0.,low)-span*.04, high+span*.12)

        for preview_ax, (label, x_min, x_max) in zip(self.preview_axes, self.visible_plot_specs()):
            self._draw_chromatogram_axis(
                preview_ax,
                x=x,
                y=y,
                y_smooth=y_smooth,
                fill_y=fill_y,
                marker_y=marker_y,
                selected_row=selected_row,
                selected_peak=selected_peak,
                title=label,
                x_min=x_min,
                x_max=x_max,
                compact=self.view_var.get() == "Все",
                normalized=True,
            )
            if self.view_var.get() == 'Пик' and selected_row is not None:
                start, end = self._selected_manual_drag_bounds()
                mask = (x >= start) & (x <= end)
                curve = next(line for line in preview_ax.lines if line.get_label() == 'Corrected')
                values = np.asarray(curve.get_ydata())[mask]
                if len(values) and np.isfinite(values).all() and np.max(values) > 0:
                    height = float(np.max(values))
                    preview_ax.set_ylim(min(0., float(np.min(values)))-height*.08, height*1.4)
        if preserve_view:
            self._restore_axes_views(axes, saved_views)
        self.draw_extra_peaks(x, y)
        self.canvas.draw_idle()
