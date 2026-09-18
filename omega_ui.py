"""Desktop workspace presentation. No chromatographic calculation lives here."""
import tkinter as tk
from tkinter import ttk
from tkinter import font as tkfont

import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from omega_core import rt_profile
from omega_version import APP_NAME

BG = '#F3F5F8'
SURFACE = '#FFFFFF'
INK = '#223047'
MUTED = '#596B80'
BLUE = '#2463C4'
LINE = '#DCE3EC'
VIEW_NAMES = ('Все', 'C16', 'C18', 'C20', 'C22', 'Пик')


def readable_status(status):
    text = str(status)
    if 'manual' in text:
        return 'Границы заданы вручную'
    if text == 'recovered_epa_local_review':
        return 'Дополнительный поиск ЭПК — проверьте вершину и границы'
    if 'not_found' in text:
        return 'Пик не найден — проверьте этот участок'
    if 'unresolved' in text:
        return 'Неоднозначное выделение — нужна проверка'
    if 'estimated' in text or 'fit' in text or 'recovered' in text:
        return 'Модельная оценка — нужна проверка'
    if text.startswith('matched'):
        return 'Пик найден автоматически'
    return 'Проверьте выделение пика'


class WorkspaceUI:
    def build_workspace(self, toolbar_class):
        style = ttk.Style(self.root)
        style.theme_use('clam')
        self.root.configure(background=BG)
        tkfont.nametofont('TkDefaultFont').configure(family='Segoe UI', size=10)
        style.configure('.', font=('Segoe UI', 9), background=BG, foreground=INK)
        style.configure('TFrame', background=BG)
        style.configure('Card.TFrame', background=SURFACE)
        style.configure('TLabel', background=BG, foreground=INK)
        style.configure('Card.TLabel', background=SURFACE)
        style.configure('Muted.TLabel', foreground=MUTED)
        style.configure('CardMuted.TLabel', background=SURFACE, foreground=MUTED, font=('Segoe UI', 9))
        style.configure('Title.TLabel', font=('Segoe UI Semibold', 15))
        style.configure('Section.TLabel', background=SURFACE, font=('Segoe UI Semibold', 11))
        style.configure('Value.TLabel', background=SURFACE, font=('Segoe UI Semibold', 24))
        style.configure('TButton', padding=(8, 5), borderwidth=0, background='#E8EDF4', foreground=INK)
        style.map('TButton', background=[('active', '#DAE4F2'), ('disabled', '#EDF0F4')],
                  foreground=[('disabled', '#8B96A5')])
        style.configure('Accent.TButton', background=BLUE, foreground='white')
        style.map('Accent.TButton', background=[('active', '#184E9C'), ('disabled', '#CAD4E3')],
                  foreground=[('disabled', '#64748B')])
        style.configure('Quiet.TButton', background=SURFACE, foreground=BLUE, padding=(4, 3))
        for name, color in [('Good', '#25714A'), ('Check', '#976014'), ('Stop', '#B53C45'), ('Pending', MUTED)]:
            tint = {'Good':'#E7F4ED', 'Check':'#FFF3D8', 'Stop':'#FCE4E7', 'Pending':'#EDF1F6'}[name]
            style.configure(f'{name}.TButton', background=tint, foreground=color, anchor='w', padding=(7, 5))
        style.configure('Tab.TRadiobutton', padding=(10, 5), background=SURFACE)
        style.layout('Tab.TRadiobutton', [('Radiobutton.padding', {'children': [('Radiobutton.label', {'sticky': 'nswe'})], 'sticky': 'nswe'})])
        style.map('Tab.TRadiobutton', background=[('selected', '#E5EEFC'), ('active', '#F0F4FA')],
                  foreground=[('selected', BLUE)])
        style.configure('TEntry', fieldbackground=SURFACE, padding=5)
        style.configure('TCombobox', padding=4, fieldbackground=SURFACE)
        style.map('TCombobox', fieldbackground=[('readonly', SURFACE)])
        style.configure('Treeview', background=SURFACE, fieldbackground=SURFACE, foreground=INK,
                        rowheight=26, borderwidth=0, font=('Segoe UI', 9))
        style.configure('Treeview.Heading', background='#F0F3F8', foreground=MUTED,
                        font=('Segoe UI Semibold', 9), padding=(5, 7), borderwidth=0)
        style.map('Treeview', background=[('selected', '#E0ECFF')], foreground=[('selected', '#174B91')])
        style.configure('TPanedwindow', background=BG)
        style.configure('Sash', sashthickness=7)
        self.view_var = tk.StringVar(value='Все')
        self.sample_search_var = tk.StringVar()
        self.review_only_var = tk.BooleanVar(value=False)
        self.review_filter_var = tk.StringVar(value='Все пробы')
        self.full_overview_var = tk.BooleanVar(value=False)
        self.sample_count_var = tk.StringVar(value='Нет загруженных проб')
        self.selected_peak_var = tk.StringVar(value='Выберите кислоту в таблице')
        self.peak_detail_var = tk.StringVar(value='Границы и сведения о пике появятся здесь.')
        self.judge_summary_var = tk.StringVar(value='Откройте CSV для расчёта')
        self.show_details_var = tk.BooleanVar(value=False)
        self.sample_title_var = tk.StringVar(value='Хроматограмма')

        header = ttk.Frame(self.root, padding=(16, 12, 16, 6))
        header.pack(fill='x')
        ttk.Label(header, text=APP_NAME, style='Title.TLabel').pack(side='left', padx=(0, 22))
        ttk.Button(header, text='Новый расчёт', style='Accent.TButton', command=self.new_calculation).pack(side='left')
        self.batch_results_button = ttk.Button(header, text='Результаты серии', command=self.open_batch_results_window)
        self.batch_results_button.pack(side='left', padx=8)
        self.tutorial_button = ttk.Button(header, text='Как работать', command=self.start_tutorial)
        self.tutorial_button.pack(side='left', padx=(0,8))
        ttk.Button(header, text='Настроить профили', command=self.open_profile_manager).pack(side='right')
        self.profile_combo = ttk.Combobox(header, textvariable=self.profile_var, state='readonly', width=17)
        self.profile_combo.pack(side='right', padx=8)
        self.profile_combo.bind('<<ComboboxSelected>>', self.on_profile_combo_selected)
        ttk.Label(header, text='Прибор', style='Muted.TLabel').pack(side='right')
        self._refresh_profile_combo()
        ttk.Label(self.root, textvariable=self.file_var, style='Muted.TLabel', padding=(16, 0, 16, 10)).pack(fill='x')

        # Reserve the footer before the expandable workspace, even on small screens.
        ttk.Label(self.root, textvariable=self.status_var, style='Muted.TLabel', padding=(16, 5)).pack(side='bottom', fill='x')
        self.main_pane = ttk.Panedwindow(self.root, orient='horizontal')
        self.main_pane.pack(fill='both', expand=True, padx=12, pady=(0, 4))
        batch_panel = ttk.Frame(self.main_pane, style='Card.TFrame', padding=10, width=215)
        plot_panel = ttk.Frame(self.main_pane, style='Card.TFrame', padding=(4, 10), width=520)
        inspector = ttk.Frame(self.main_pane, style='Card.TFrame', padding=12, width=355)
        self.main_pane.add(batch_panel, weight=0)
        self.main_pane.add(plot_panel, weight=1)
        self.main_pane.add(inspector, weight=0)
        self.inspector = inspector

        ttk.Label(batch_panel, text='Пробы', style='Section.TLabel').pack(anchor='w', pady=(0, 10))
        ttk.Label(batch_panel, text='Номер или имя пробы', style='CardMuted.TLabel').pack(anchor='w')
        self.search_entry = ttk.Entry(batch_panel, textvariable=self.sample_search_var, width=18)
        self.search_entry.pack(fill='x', pady=(4, 8))
        filters = ttk.Combobox(batch_panel, textvariable=self.review_filter_var, state='readonly', width=18,
                              values=('Все пробы','Требуют проверки','Требуют обязательной проверки'))
        filters.pack(fill='x', pady=(0, 8))
        filters.bind('<<ComboboxSelected>>', lambda e: self.populate_main_batch_tree())
        self.filter_description = ttk.Label(batch_panel, text='', style='CardMuted.TLabel', wraplength=190)
        self.filter_description.pack(fill='x')
        nav = ttk.Frame(batch_panel, style='Card.TFrame')
        nav.pack(fill='x', pady=(0, 8))
        self.prev_button = ttk.Button(nav, text='Назад', command=self.prev_batch, width=7)
        self.prev_button.pack(side='left')
        self.next_button = ttk.Button(nav, text='Далее', command=self.next_batch, width=7)
        self.next_button.pack(side='right')
        ttk.Label(batch_panel, textvariable=self.batch_var, style='CardMuted.TLabel').pack(anchor='w', pady=(0, 8))
        ttk.Label(batch_panel, textvariable=self.sample_count_var, style='CardMuted.TLabel').pack(side='bottom', anchor='w', pady=(8, 0))
        batch_host = ttk.Frame(batch_panel, style='Card.TFrame')
        batch_host.pack(fill='both', expand=True)
        self.batch_tree = ttk.Treeview(batch_host, columns=('sample_name','omega_value','confidence'),
                                      show='headings', selectmode='browse', height=8)
        for name, label, width in [('sample_name','Проба',45), ('omega_value','%',42), ('confidence','Статус',85)]:
            self.batch_tree.heading(name, text=label)
            self.batch_tree.column(name, width=width, minwidth=width, stretch=name=='sample_name', anchor='e' if name=='omega_value' else 'w')
        self.batch_tree.pack(side='left', fill='both', expand=True)
        scrollbar = ttk.Scrollbar(batch_host, orient='vertical', command=self.batch_tree.yview)
        scrollbar.pack(side='right', fill='y')
        self.batch_tree.configure(yscrollcommand=scrollbar.set)
        self.batch_tree.bind('<<TreeviewSelect>>', self.handle_batch_tree_selection)
        self._configure_quality_tags(self.batch_tree)
        self.sample_search_var.trace_add('write', lambda *_: self.populate_main_batch_tree())

        ttk.Label(plot_panel, textvariable=self.sample_title_var, style='Section.TLabel', padding=(10, 0)).pack(anchor='w')
        views = ttk.Frame(plot_panel, style='Card.TFrame', padding=(8, 8, 0, 4))
        views.pack(fill='x')
        for name in VIEW_NAMES:
            ttk.Radiobutton(views, text=name, value=name, variable=self.view_var, style='Tab.TRadiobutton',
                            command=self.change_plot_view).pack(side='left', padx=(0, 3))
        ttk.Checkbutton(plot_panel, text='Показать всю запись, включая поздние пики', variable=self.full_overview_var,
                        command=lambda: self.update_plot() if self.df_processed is not None else None).pack(anchor='w', padx=8)
        toolbar_host = tk.Frame(plot_panel, background=SURFACE)
        toolbar_host.pack(side='bottom', fill='x')
        self.figure = Figure(figsize=(5, 4), dpi=100, facecolor=SURFACE)
        self.preview_specs = rt_profile.preview_windows(self.reference_targets)
        self._make_plot_axes()
        self.canvas = FigureCanvasTkAgg(self.figure, master=plot_panel)
        self.canvas.get_tk_widget().configure(width=300, height=260)
        self.canvas.get_tk_widget().pack(fill='both', expand=True)
        for event, callback in [('button_press_event', self.handle_manual_boundary_press),
                                ('motion_notify_event', self.handle_manual_boundary_motion),
                                ('button_release_event', self.handle_manual_boundary_release)]:
            self.canvas.mpl_connect(event, callback)
        self.canvas.mpl_connect('key_press_event', lambda event: self.cancel_manual_pick() if event.key=='escape' else None)
        self.plot_toolbar = toolbar_class(self.canvas, toolbar_host, pack_toolbar=False)
        self.plot_toolbar.configure(background=SURFACE)
        self.plot_toolbar.update()
        self.plot_toolbar.pack(side='left', fill='x', expand=True)

        ttk.Label(inspector, text='ОМЕГА-3 ИНДЕКС', style='CardMuted.TLabel').pack(anchor='w')
        ttk.Label(inspector, textvariable=self.omega_var, style='Value.TLabel').pack(anchor='w', pady=(0, 2))
        self.judge_label = ttk.Button(inspector, textvariable=self.judge_summary_var,
                                     style='Pending.TButton', command=self.show_confidence_details)
        self.judge_label.pack(fill='x')
        self.review_buttons = ttk.Frame(inspector, style='Card.TFrame')
        self.review_buttons.pack(fill='x', pady=(4, 0))
        self.confidence_button = self.judge_label
        self.confidence_button.state(['disabled'])
        ttk.Separator(inspector).pack(fill='x', pady=(0, 8))
        # Bottom packing guarantees access to manual controls when the table shrinks.
        manual = ttk.Frame(inspector, style='Card.TFrame')
        manual.pack(side='bottom', fill='x')
        self.manual_panel = manual
        ttk.Separator(manual).pack(fill='x', pady=(8, 8))
        ttk.Label(manual, text='Ручная интеграция', style='Section.TLabel').pack(anchor='w')
        self.peak_name_label = ttk.Label(manual, textvariable=self.selected_peak_var, style='Card.TLabel', wraplength=310)
        self.peak_name_label.pack(fill='x', pady=(4, 2))
        self.peak_detail_label = ttk.Label(manual, textvariable=self.peak_detail_var, style='CardMuted.TLabel', wraplength=310)
        self.peak_detail_label.pack(fill='x', pady=(0, 4))
        bounds = ttk.Frame(manual, style='Card.TFrame')
        bounds.pack(fill='x')
        bounds.columnconfigure((0, 1), weight=1, uniform='bounds')
        for col, label, var in [(0, 'Начало, мин', self.manual_start_var), (1, 'Конец, мин', self.manual_end_var)]:
            ttk.Label(bounds, text=label, style='CardMuted.TLabel').grid(row=0, column=col, sticky='w', padx=(0, 6))
            ttk.Entry(bounds, textvariable=var, width=12).grid(row=1, column=col, sticky='ew', padx=(0, 6) if col==0 else 0, pady=(3, 6))
        selection_actions = ttk.Frame(manual, style='Card.TFrame')
        selection_actions.pack(fill='x', pady=(0, 6))
        self.manual_pick_button = ttk.Button(selection_actions, text='Выделить границы', command=self.toggle_manual_pick)
        self.manual_pick_button.pack(side='left', fill='x', expand=True)
        self.clear_selection_button = ttk.Button(selection_actions, text='Отмена выделения', command=self.clear_peak_selection)
        self.clear_selection_button.pack(side='left', fill='x', expand=True, padx=(6, 0))
        actions = ttk.Frame(manual, style='Card.TFrame')
        actions.pack(fill='x')
        self.apply_bounds_button = ttk.Button(actions, text='Применить границы', style='Accent.TButton', command=self.apply_manual_integration)
        self.apply_bounds_button.pack(side='left', fill='x', expand=True)
        ttk.Button(actions, text='Отменить', command=self.undo_manual_edit).pack(side='right', padx=(6, 0))
        transfers = ttk.Frame(manual, style='Card.TFrame')
        transfers.pack(fill='x', pady=(6, 0))
        ttk.Button(transfers, text='Сохранить правки', style='Quiet.TButton', command=self.save_manual_edits).pack(side='left')
        ttk.Button(transfers, text='Загрузить правки', style='Quiet.TButton', command=self.load_manual_edits).pack(side='right')
        self.extra_peaks_var = tk.StringVar(value='')
        ttk.Label(manual, textvariable=self.extra_peaks_var, style='CardMuted.TLabel', wraplength=310).pack(anchor='w')

        peaks_header = ttk.Frame(inspector, style='Card.TFrame')
        peaks_header.pack(fill='x', pady=(0, 5))
        ttk.Label(peaks_header, text='Кислоты', style='Section.TLabel').pack(side='left')
        ttk.Checkbutton(peaks_header, text='Подробности', variable=self.show_details_var, command=self.toggle_peak_columns).pack(side='right')
        table = ttk.Frame(inspector, style='Card.TFrame')
        table.pack(fill='both', expand=True)
        cols = ('display_name','area','percent_area','code','expected_rt','found_rt','status')
        self.tree = ttk.Treeview(table, columns=cols, displaycolumns=('display_name','area','percent_area'),
                                show='headings', selectmode='browse', height=5)
        for key, label, width in zip(cols, ('Кислота','Площадь','%','Код','RT расч.','RT найден','Комментарий'), (150,76,52,82,78,78,240)):
            self.tree.heading(key, text=label)
            self.tree.column(key, width=width, minwidth=width, stretch=key=='display_name', anchor='w' if key in ('display_name','status','code') else 'e')
        self.tree.grid(row=0, column=0, sticky='nsew')
        scroll_y = ttk.Scrollbar(table, orient='vertical', command=self.tree.yview)
        scroll_y.grid(row=0, column=1, sticky='ns')
        scroll_x = ttk.Scrollbar(table, orient='horizontal', command=self.tree.xview)
        scroll_x.grid(row=1, column=0, sticky='ew')
        self.tree.configure(yscrollcommand=scroll_y.set, xscrollcommand=scroll_x.set)
        table.columnconfigure(0, weight=1)
        table.rowconfigure(0, weight=1)
        self.tree.bind('<<TreeviewSelect>>', self.handle_target_selection)
        self.tree.tag_configure('review', foreground='#976014')
        self.inspector.bind('<Configure>', self._resize_inspector)
        self.update_batch_navigation()
        self.root.after_idle(self._apply_responsive_pane_positions)

    def _resize_inspector(self, event):
        width = max(170, event.width-28)
        for widget in (self.peak_name_label, self.peak_detail_label):
            widget.configure(wraplength=width)
        compact = event.height < 690
        if compact != getattr(self, '_compact_workspace', None):
            self._compact_workspace = compact
            self.update_peak_details()
        code = self.selected_target_code
        if code and self.tree.exists(code):
            self.root.after_idle(lambda: self.tree.see(code) if self.tree.exists(code) else None)

    def position_workspace(self):
        self.root.update_idletasks()
        width = self.main_pane.winfo_width()
        if width > 1:
            self.main_pane.sashpos(0, 230)
            self.main_pane.sashpos(1, width-(370 if width>1200 else 350))

    def _make_plot_axes(self):
        self.figure.clear()
        self._manual_overlay_artists = {}
        mode = self.view_var.get()
        if mode == 'Все':
            grid = self.figure.add_gridspec(3, 2, height_ratios=[.9, 1.35, 1.35], hspace=.66, wspace=.25)
            self.ax = self.figure.add_subplot(grid[0, :])
            self.preview_axes = [self.figure.add_subplot(grid[1+i//2, i%2]) for i in range(4)]
        else:
            grid = self.figure.add_gridspec(2, 1, height_ratios=[1, 3.5], hspace=.44)
            self.ax = self.figure.add_subplot(grid[0])
            self.preview_axes = [self.figure.add_subplot(grid[1])]
        self.figure.subplots_adjust(left=.11, right=.97, top=.93, bottom=.12)
        for axis in (self.ax, *self.preview_axes):
            axis.set_facecolor(SURFACE)
            axis.tick_params(labelsize=8, colors=MUTED)
            for spine in axis.spines.values():
                spine.set_color(LINE)
        if self.df_processed is None:
            for axis in (self.ax, *self.preview_axes):
                axis.set_axis_off()
            self.ax.text(.5, .5, 'Откройте CSV, чтобы увидеть хроматограмму', transform=self.ax.transAxes,
                         ha='center', va='center', color=MUTED, fontsize=10)

    def change_plot_view(self):
        self.cancel_manual_pick()
        self._make_plot_axes()
        self.plot_toolbar.update()
        if self.df_processed is not None:
            self.update_plot()
        else:
            self.canvas.draw_idle()

    def visible_plot_specs(self):
        mode = self.view_var.get()
        if mode == 'Все':
            return [(f'{group} · {label} мин', a, b) for group, (label, a, b) in zip(VIEW_NAMES[1:5], self.preview_specs)]
        if mode in VIEW_NAMES[1:5]:
            label, a, b = self.preview_specs[VIEW_NAMES.index(mode)-1]
            return [(f'{mode} · {label} мин', a, b)]
        rows = self.matched_targets_df.loc[self.matched_targets_df.code == self.selected_target_code]
        if not rows.empty:
            row = rows.iloc[0]
            center = row.found_rt if np.isfinite(row.found_rt) else row.get('corrected_target_rt', row.expected_rt)
            a, b = row.integration_start_x, row.integration_end_x
            if not np.isfinite([a, b]).all():
                a, b = center-.025, center+.025
            margin = max(.025, (b-a)*.65)
            return [(f'{row.code} · масштаб выбранного пика', a-margin, b+margin)]
        label, a, b = self.preview_specs[0]
        return [('Выберите кислоту в таблице', a, b)]

    def clear_peak_selection(self):
        self.cancel_manual_pick()
        if self._manual_drag_after_id is not None:
            self.root.after_cancel(self._manual_drag_after_id)
            self._manual_drag_after_id = None
        self._manual_drag_active_boundary = None
        self._manual_drag_pending_bounds = None
        self._manual_drag_axis = None
        self.selected_target_code = None
        self.tree.selection_remove(*self.tree.selection())
        self.tree.focus('')
        self.load_selected_integration_bounds(silent=True)
        self.update_peak_details()
        if self.view_var.get() == 'Пик':
            self.view_var.set('Все')
            self.change_plot_view()
        elif self.df_processed is not None:
            self.update_plot(preserve_view=True)
        self.status_var.set('Выделение кислоты снято. Применённые правки сохранены.')

    def focus_peak(self, code):
        if not self.tree.exists(code):
            return
        already_selected = code == self.selected_target_code
        self.tree.selection_set(code)
        self.tree.focus(code)
        self.tree.see(code)
        self.handle_target_selection()
        if already_selected:
            self.view_var.set(next((g for g in VIEW_NAMES[1:5] if code.startswith(g+':')), 'Пик'))
            self.change_plot_view()

    def update_peak_details(self):
        if self.matched_targets_df.empty or not self.selected_target_code:
            self.selected_peak_var.set('Выберите кислоту в таблице')
            self.peak_detail_var.set('Щёлкните по кислоте, чтобы открыть её участок.')
            return
        rows = self.matched_targets_df.loc[self.matched_targets_df.code == self.selected_target_code]
        if rows.empty:
            return
        row = rows.iloc[0]
        self.selected_peak_var.set(f"{row.display_name} · {row.code}")
        rt = f'RT {row.found_rt:.4f} мин · ' if np.isfinite(row.found_rt) else ''
        if getattr(self, '_compact_workspace', False):
            status = 'Вручную' if 'manual' in str(row.status) else 'Проверить' if 'recovered' in str(row.status) or 'unresolved' in str(row.status) else 'Автоматически' if str(row.status).startswith('matched') else 'Не найден'
            self.peak_detail_var.set(rt + status)
        else:
            self.peak_detail_var.set(rt + readable_status(row.status))

    def update_review_panel(self):
        confidence = self.current_confidence or {}
        text, tag = self._operator_quality(confidence)
        titles = {'quality_good': 'Проверки пройдены', 'quality_check': 'Требуется проверка',
                  'quality_stop': 'Стоп — проверьте пики', 'quality_pending': 'Ожидание расчёта'}
        styles = {'quality_good': 'Good', 'quality_check': 'Check', 'quality_stop': 'Stop', 'quality_pending': 'Pending'}
        self.judge_summary_var.set(titles[tag] + ' · подробнее')
        self.judge_label.configure(style=f'{styles[tag]}.TButton')
        for child in self.review_buttons.winfo_children():
            child.destroy()
        codes = confidence.get('high_error_risk', {}).get('peak_codes', [])
        aliases = {'C20:5': 'ЭПК', 'C22:6': 'ДГК'}
        for index, code in enumerate(codes):
            ttk.Button(self.review_buttons, text=aliases.get(code, code), style='Quiet.TButton', width=7,
                       command=lambda c=code: self.focus_peak(c)).grid(row=index//4, column=index%4, sticky='w', padx=(0, 4))
        self.update_peak_details()

    def toggle_peak_columns(self):
        self.tree.configure(displaycolumns=('#all' if self.show_details_var.get() else ('display_name','area','percent_area')))
