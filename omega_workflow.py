"""Profile selection, explicit peak insertion and optional exit review."""
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import numpy as np
from omega_core import instrument_profiles, manual_edits, integration


class WorkflowUI:
    def setup_workflow(self):
        self._insert_pick = None
        self._new_profile_dialog = None
        self.quiz_window = None
        self.tutorial_window = None
        self.root.protocol('WM_DELETE_WINDOW', self.request_close)
        self.canvas.mpl_connect('button_press_event', self.context_peak_menu)
        menu = tk.Menu(self.root)
        training = tk.Menu(menu, tearoff=False)
        training.add_command(label='Квиз по пикам', command=lambda:self.start_peak_quiz(False))
        training.add_command(label='Импортировать ответы квиза…', command=self.import_quiz_answers)
        training.add_command(label='Показать папку ответов', command=self.show_quiz_folder)
        menu.add_cascade(label='Обучение', menu=training)
        help_menu = tk.Menu(menu, tearoff=False)
        help_menu.add_command(label='Как работать с Omega', command=self.start_tutorial)
        menu.add_cascade(label='Помощь', menu=help_menu)
        self.root.configure(menu=menu)

    def start_application(self):
        from omega_tutorial import should_show_tutorial, mark_tutorial_seen
        if should_show_tutorial():
            def finish_onboarding():
                mark_tutorial_seen()
                self.new_calculation()
            self.start_tutorial(on_done=finish_onboarding)
        else:
            self.new_calculation()

    def start_tutorial(self, on_done=None):
        if self._background_future is not None or self._preload_after_id is not None:
            self.status_var.set('Обучалку можно открыть после завершения расчёта.')
            return
        if self.tutorial_window is not None and self.tutorial_window.window.winfo_exists():
            self.tutorial_window.window.lift()
            return
        from omega_tutorial import OperatorTutorial
        self.tutorial_window = OperatorTutorial(self, on_done=on_done)

    def preserve_manual_work(self):
        if not any(b.get('manual_overrides') or b.get('manual_extra_peaks') for b in self.loaded_batches):
            return True
        answer=messagebox.askyesnocancel('Ручные правки','Сохранить ручные правки перед продолжением?',parent=self.root)
        if answer is None:
            return False
        if answer:
            return bool(self.save_manual_edits())
        return True

    def new_calculation(self):
        if self._background_future is not None or self._preload_after_id is not None:
            return
        if self._new_profile_dialog is not None and self._new_profile_dialog.winfo_exists():
            self._new_profile_dialog.lift()
            return
        if not self.preserve_manual_work():
            return
        window=tk.Toplevel(self.root)
        self._new_profile_dialog=window
        window.title('Новый расчёт — профиль прибора')
        window.geometry('540x280')
        window.transient(self.root)
        form=ttk.Frame(window,padding=22)
        form.pack(fill='both',expand=True)
        ttk.Label(form,text='Выберите профиль для нового CSV',font=('Segoe UI',13,'bold')).pack(anchor='w',pady=(0,12))
        self.new_profile_var=tk.StringVar(value=self.active_instrument_profile['name'])
        combo=ttk.Combobox(form,textvariable=self.new_profile_var,state='readonly',width=40,
                           values=[p['name'] for p in self.profile_store['profiles']])
        combo.pack(fill='x',pady=6)
        ttk.Label(form,text='Профиль будет применён только после выбора и чтения CSV.',wraplength=470).pack(anchor='w',pady=8)
        def configure():
            self.open_profile_manager()
            window.wait_window(self.profile_window)
            combo.configure(values=[p['name'] for p in self.profile_store['profiles']])
            window.lift()
        ttk.Button(form,text='Настроить или создать профиль',command=configure).pack(anchor='w')
        def close():
            window.destroy()
            self._new_profile_dialog=None
        def proceed():
            profile=next((p for p in self.profile_store['profiles'] if p['name']==self.new_profile_var.get()),None)
            if profile is None:
                return
            close()
            self.open_file(profile_id=profile['id'])
        bar=ttk.Frame(form)
        bar.pack(fill='x',pady=(18,0))
        ttk.Button(bar,text='Отмена',command=close).pack(side='left')
        ttk.Button(bar,text='Выбрать CSV и рассчитать',style='Accent.TButton',command=proceed).pack(side='right')
        window.protocol('WM_DELETE_WINDOW',close)
        window.bind('<Return>',lambda e:proceed())

    def finish_new_calculation(self, path, batches, profile_id):
        if not batches:
            raise ValueError('В CSV нет проб для анализа.')
        if profile_id is not None:
            if self._profile_by_id(profile_id) is None:
                raise ValueError('Выбранный профиль больше не существует.')
            self.profile_store['active_profile_id']=profile_id
            self.profile_store=instrument_profiles.save_store(self.profile_store,self.base_reference_targets)
            self.active_instrument_profile=instrument_profiles.active_profile(self.profile_store)
            self.reference_targets=instrument_profiles.apply_profile_to_targets(self.base_reference_targets,self.active_instrument_profile)
            self._refresh_profile_combo()
        self._finish_file_load(path,batches)

    def context_peak_menu(self,event):
        if event.button!=3 or event.inaxes not in [self.ax,*self.preview_axes] or self.df_processed is None:
            return
        menu=tk.Menu(self.root,tearoff=False)
        menu.add_command(label='Добавить пик — указать две границы',command=self.begin_peak_insertion)
        for peak in self.loaded_batches[self.current_batch_index].get('manual_extra_peaks',[]):
            if event.xdata is not None and peak['start']<=event.xdata<=peak['end']:
                menu.add_command(label=f"Назначить кислоту: {peak['name']}",
                                 command=lambda p=peak:self.insert_peak_dialog(p['start'],p['end']))
        if self.selected_target_code:
            menu.add_command(label=f'Изменить границы {self.selected_target_code}',command=self.toggle_manual_pick)
        menu.tk_popup(self.root.winfo_pointerx(),self.root.winfo_pointery())
        menu.grab_release()

    def begin_peak_insertion(self):
        self.cancel_manual_pick()
        if self.plot_toolbar.mode=='pan/zoom':self.plot_toolbar.pan()
        elif self.plot_toolbar.mode=='zoom rect':self.plot_toolbar.zoom()
        self._insert_pick=dict(first=None,batch=self.current_batch_index)
        self.status_var.set('Новый пик: щёлкните левую и правую границу. Esc — отмена.')

    def insertion_click(self,event):
        if self._insert_pick is None:return False
        if event.button!=1 or event.inaxes not in [self.ax,*self.preview_axes] or event.xdata is None:return True
        if self._insert_pick['batch']!=self.current_batch_index:
            self.cancel_manual_pick();return True
        x=self.df_processed.x_corrected.to_numpy()
        value=float(event.xdata)
        if not x[0]<=value<=x[-1]:return True
        if self._insert_pick['first'] is None:
            self._insert_pick['first']=value
            self._manual_pick_artist=event.inaxes.axvline(value,color='#7C3AED',linestyle='--')
            self.canvas.draw_idle()
        else:
            a,b=sorted((value,self._insert_pick['first']))
            self.cancel_manual_pick()
            if a<b:self.insert_peak_dialog(a,b)
        return True

    def insert_peak_dialog(self,start,end):
        window=tk.Toplevel(self.root);window.title('Добавить пик');window.transient(self.root)
        form=ttk.Frame(window,padding=18);form.pack(fill='both',expand=True)
        ttk.Label(form,text='Назначение пика',font=('Segoe UI',12,'bold')).pack(anchor='w')
        choices={'Не назначен — не включать в Омега-3':None}
        choices.update({f'{r.code} · {r.display_name}':r.code for r in self.matched_targets_df.itertuples()})
        target=tk.StringVar(value=next(iter(choices)))
        ttk.Combobox(form,textvariable=target,values=list(choices),state='readonly',width=55).pack(fill='x',pady=10)
        a,b=tk.StringVar(value=f'{start:.8f}'),tk.StringVar(value=f'{end:.8f}')
        for label,var in [('Начало, мин',a),('Конец, мин',b)]:
            ttk.Label(form,text=label).pack(anchor='w');ttk.Entry(form,textvariable=var).pack(fill='x',pady=4)
        note=ttk.Label(form,text='Неназначенный пик сохраняется отдельно от расчёта. Чтобы учесть кислоту, выберите её из списка.',wraplength=440)
        note.pack(fill='x',pady=8)
        def save():
            try:
                lo,hi=float(a.get().replace(',','.')),float(b.get().replace(',','.'))
                code=choices[target.get()]
                batch=self.loaded_batches[self.current_batch_index]
                if code:
                    row=self.matched_targets_df.loc[self.matched_targets_df.code==code].iloc[0]
                    if np.isfinite(row.found_rt) and not messagebox.askyesno('Заменить выделение?',f'{code} уже найден. Применить к нему новые границы?',parent=window):return
                    manual_edits.apply_edit(batch,code,lo,hi,source='insert')
                    self.selected_target_code=code
                else:
                    manual_edits.add_unassigned_peak(batch,lo,hi,f"Пик {len(batch.get('manual_extra_peaks',[]))+1}")
                self._refresh_after_manual_edit()
                window.destroy()
            except (ValueError,KeyError) as exc:
                messagebox.showerror('Пик не добавлен',str(exc),parent=window)
        ttk.Button(form,text='Добавить пик',style='Accent.TButton',command=save).pack(fill='x',pady=8)

    def draw_extra_peaks(self,x,y):
        if not self.loaded_batches:return
        extras=self.loaded_batches[self.current_batch_index].get('manual_extra_peaks',[])
        self.extra_peaks_var.set(f'Неназначенных пиков: {len(extras)} · вне расчёта' if extras else '')
        for axis in [self.ax,*self.preview_axes]:
            curve=next((line for line in axis.lines if line.get_label()=='Corrected'),None)
            if curve is None:continue
            drawn=np.asarray(curve.get_ydata())
            for peak in extras:
                if peak['end']<axis.get_xlim()[0] or peak['start']>axis.get_xlim()[1]:continue
                fx,fy=integration.interval_points(x,drawn,peak['start'],peak['end'])
                axis.fill_between(fx,0,np.maximum(fy,0),color='#7C3AED',alpha=.22)
                axis.axvline(peak['start'],color='#7C3AED',ls='--',lw=.8)
                axis.axvline(peak['end'],color='#7C3AED',ls='--',lw=.8)
                axis.annotate(peak['name']+' (не назначен)',xy=(peak['apex'],np.interp(peak['apex'],x,drawn)),
                              xytext=(4,12),textcoords='offset points',color='#6D28D9',fontsize=8)

    def show_quiz_folder(self):
        from omega_review.quiz_data import outbox
        messagebox.showinfo('Файлы ответов',str(outbox()),parent=self.root)

    def import_quiz_answers(self):
        from omega_review.quiz_data import import_package
        path=filedialog.askopenfilename(parent=self.root,filetypes=[('Ответы квиза Omega','*.omega-quiz')])
        if not path:return
        try:
            count=import_package(path)
            messagebox.showinfo('Импорт завершён',f'Добавлено ответов: {count}. Обучение запускается отдельно.',parent=self.root)
        except Exception as exc:
            messagebox.showerror('Не удалось импортировать',str(exc),parent=self.root)

    def request_close(self):
        if self._background_future is not None or self._preload_after_id is not None:
            messagebox.showinfo('Идёт расчёт','Дождитесь окончания расчёта перед выходом.',parent=self.root);return
        if self.quiz_window is not None and self.quiz_window.window.winfo_exists():
            self.quiz_window.window.lift();return
        if not self.preserve_manual_work():return
        from omega_review import data
        available=bool(self.loaded_batches) or any((data.workspace()/'examples').glob('*/snapshot.json'))
        if not available:self.root.destroy();return
        answer=messagebox.askyesnocancel('Квиз по пикам','Перед выходом оценить до 6 пиков для обучения судьи?\nДа — пройти квиз; Нет — выйти; Отмена — продолжить работу.',parent=self.root)
        if answer is None:return
        if answer:self.start_peak_quiz(True)
        else:self.root.destroy()

    def start_peak_quiz(self,close_after=False):
        from omega_review.quiz_data import choose_questions
        from omega_review.quiz_ui import PeakQuiz
        if self.quiz_window is not None and self.quiz_window.window.winfo_exists():
            self.quiz_window.window.lift();return
        try:
            questions=choose_questions(self)
            if not questions:
                messagebox.showinfo('Квиз','Сначала загрузите и рассчитайте CSV.',parent=self.root)
                if close_after:self.root.destroy()
                return
            self.quiz_window=PeakQuiz(self,questions,self.root.destroy if close_after else lambda:None)
        except Exception as exc:
            messagebox.showerror('Квиз недоступен',str(exc),parent=self.root)
