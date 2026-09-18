"""Optional six-peak exit quiz; records only answers explicitly submitted."""
import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.ticker import ScalarFormatter, MaxNLocator
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
from . import data, quiz_data


class PeakQuiz:
    def __init__(self, app, questions, on_done, root=None):
        self.app, self.questions, self.on_done, self.storage = app, questions, on_done, root
        self.index, self.answers = 0, []
        self.window = tk.Toplevel(app.root)
        self.window.title('Квиз по пикам')
        self.window.geometry('980x680')
        self.window.minsize(820,600)
        self.window.transient(app.root)
        self.window.protocol('WM_DELETE_WINDOW',self.finish)
        self.title = tk.StringVar()
        ttk.Label(self.window,textvariable=self.title,font=('Segoe UI',14,'bold'),padding=12).pack(anchor='w')
        ttk.Label(self.window,text='Оцените границы и назначение пика. Ответов модели и прежних оценок здесь нет.',padding=(12,0)).pack(anchor='w')
        actions = ttk.Frame(self.window,padding=12)
        actions.pack(side='bottom',fill='x')
        ttk.Button(actions,text='Пропустить',command=self.skip).pack(side='left')
        ttk.Button(actions,text='Завершить сейчас',command=self.finish).pack(side='left',padx=8)
        self.next_button = ttk.Button(actions,text='Сохранить ответ и далее',style='Accent.TButton',command=self.answer)
        self.next_button.pack(side='right')
        self.location = ttk.Label(self.window,text=f'Файлы ответов: {quiz_data.outbox(root)}',wraplength=900,padding=(12,6))
        self.location.pack(side='bottom',fill='x')
        body = ttk.Frame(self.window,padding=12)
        body.pack(fill='both',expand=True)
        form = ttk.Frame(body,width=270)
        form.pack(side='right',fill='y',padx=(14,0))
        self.verdict = tk.StringVar(value='Верно')
        ttk.Label(form,text='Разметка пика').pack(anchor='w')
        ttk.Combobox(form,textvariable=self.verdict,values=('Верно','Есть замечания','Не уверен'),state='readonly',width=26).pack(fill='x',pady=5)
        self.issues = {}
        for key,label in data.ISSUES.items():
            var=tk.BooleanVar(value=False)
            self.issues[key]=var
            ttk.Checkbutton(form,text=label,variable=var,command=lambda:self.verdict.set('Есть замечания')).pack(anchor='w',pady=2)
        self.serious=tk.BooleanVar(value=False)
        ttk.Checkbutton(form,text='Серьёзная ошибка',variable=self.serious,command=lambda:self.verdict.set('Есть замечания')).pack(anchor='w',pady=5)
        ttk.Label(form,text='Оценка этого пика (необязательно)').pack(anchor='w',pady=(8,0))
        self.grade=tk.StringVar(value='Без оценки')
        ttk.Combobox(form,textvariable=self.grade,values=('Без оценки','5','4','3','2','1','0'),state='readonly',width=15).pack(anchor='w',pady=4)
        ttk.Label(form,text='Комментарий').pack(anchor='w')
        self.note=tk.Text(form,height=3,width=28,wrap='word')
        self.note.pack(fill='x')
        plot=ttk.Frame(body)
        plot.pack(fill='both',expand=True)
        self.figure=Figure(figsize=(5,4),dpi=100,layout='constrained')
        self.axis=self.figure.add_subplot()
        self.canvas=FigureCanvasTkAgg(self.figure,master=plot)
        self.toolbar=NavigationToolbar2Tk(self.canvas,plot,pack_toolbar=False)
        self.toolbar.pack(side='bottom',fill='x')
        self.canvas.get_tk_widget().pack(fill='both',expand=True)
        self.show_question()

    def show_question(self):
        if self.index >= len(self.questions):
            self.next_button.configure(text='Сохранить и завершить',command=self.finish)
            self.finish()
            return
        q=self.questions[self.index]
        row=next(r for r in q['snapshot']['metadata']['rows'] if r['code']==q['code'])
        name=row.get('display_name') or q['code']
        self.title.set(f"Пик {self.index+1} из {len(self.questions)} · {name} · {q['code']}")
        self.verdict.set('Верно');self.grade.set('Без оценки');self.serious.set(False)
        for var in self.issues.values():var.set(False)
        self.note.delete('1.0','end')
        signal=q['snapshot']['signal'];x=signal[:,0];y=signal[:,1]-signal[:,2]
        center=data.finite(row.get('found_rt'),data.finite(row.get('corrected_target_rt'),row['expected_rt']))
        a=data.finite(row.get('integration_start_x'),center-.02)
        b=data.finite(row.get('integration_end_x'),center+.02)
        # Include the complete integration intervals of the nearest local
        # neighbours, rather than showing only their clipped flanks.
        lo,hi=min(center-.08,a-.03),max(center+.08,b+.03)
        neighbours=[]
        for other in q['snapshot']['metadata']['rows']:
            rt=data.finite(other.get('found_rt'),np.nan)
            if other['code']!=q['code'] and np.isfinite(rt) and 0<abs(rt-center)<=.3:
                neighbours.append((rt,other))
        left=[item for item in neighbours if item[0]<center]
        right=[item for item in neighbours if item[0]>center]
        adjacent=([max(left,key=lambda item:item[0])] if left else [])
        adjacent+=([min(right,key=lambda item:item[0])] if right else [])
        for rt,other in adjacent:
            lo=min(lo,rt-.03,data.finite(other.get('integration_start_x'),rt)-.03)
            hi=max(hi,rt+.03,data.finite(other.get('integration_end_x'),rt)+.03)
        # Unassigned neighbours also need to fit. Smoothing is used only to
        # choose the view; the plotted and reviewed signal stays unsmoothed.
        local=(x>=center-.3)&(x<=center+.3)
        lx,ly=x[local],y[local]
        if len(lx)>3:
            smooth=gaussian_filter1d(ly,max(.5,.001/float(np.median(np.diff(lx)))))
            peaks,_=find_peaks(smooth,prominence=max(float(np.ptp(smooth))*.015,1e-8))
            left=peaks[lx[peaks]<min(a,center)]
            right=peaks[lx[peaks]>max(b,center)]
            if len(left):
                apex=int(left[-1])
                previous=int(left[-2]) if len(left)>1 else 0
                valley=previous+int(np.argmin(smooth[previous:apex+1]))
                lo=min(lo,float(lx[valley]))
            if len(right):
                apex=int(right[0])
                following=int(right[1]) if len(right)>1 else len(lx)-1
                valley=apex+int(np.argmin(smooth[apex:following+1]))
                hi=max(hi,float(lx[valley]))
        lo,hi=max(float(x[0]),lo),min(float(x[-1]),hi)
        mask=(x>=lo)&(x<=hi)
        self.axis.clear()
        self.axis.plot(x[mask],y[mask],color='#172333',lw=1)
        if row.get('integration_start_x') is not None and row.get('integration_end_x') is not None:
            from omega_core.integration import interval_points
            fx,fy=interval_points(x,y,a,b)
            self.axis.fill_between(fx,0,np.maximum(fy,0),color='#2463C4',alpha=.25)
            self.axis.axvline(a,color='#2463C4');self.axis.axvline(b,color='#2463C4')
        if mask.any():
            bottom=min(0.,float(np.min(y[mask])))
            top=max(0.,float(np.max(y[mask])))
            padding=max(top-bottom,1e-8)
            self.axis.set_ylim(bottom-.06*padding,top+.18*padding)
        if x[0]<=center<=x[-1]:
            height=float(np.interp(center,x,y))
            self.axis.plot(center,height,'o',color='#2463C4',ms=5)
            self.axis.annotate(q['code'],xy=(center,height),xytext=(0,14),
                textcoords='offset points',ha='center',color='#2463C4',fontweight='bold',
                bbox=dict(facecolor='white',edgecolor='none',alpha=.9,pad=2))
        self.axis.set_xlim(lo,hi)
        self.axis.set_xlabel('Время, мин');self.axis.set_ylabel('Сигнал')
        self.axis.set_title(f"Оцениваем: {name} · {q['code']}\nНесглаженный сигнал · время {center:.4f} мин",fontsize=10)
        self.axis.grid(alpha=.2)
        self.axis.xaxis.set_major_locator(MaxNLocator(6))
        self.axis.xaxis.set_major_formatter(ScalarFormatter(useOffset=False))
        self.toolbar.update()
        self.canvas.draw_idle()

    def answer(self):
        q=self.questions[self.index]
        verdict={'Верно':'good','Есть замечания':'bad','Не уверен':'unsure'}[self.verdict.get()]
        issues=[key for key,var in self.issues.items() if var.get()]
        if verdict=='bad' and not issues:issues=['other']
        label=dict(verdict=verdict,issues=issues,note=self.note.get('1.0','end').strip(),serious=self.serious.get())
        try:
            data.validate_review(q['snapshot'],None,{q['code']:label})
        except ValueError as exc:
            messagebox.showerror('Проверьте ответ',str(exc),parent=self.window)
            return
        grade=None if self.grade.get()=='Без оценки' else int(self.grade.get())
        self.answers.append({**q,'label':label,'peak_rating':grade})
        self.index+=1
        self.show_question()

    def skip(self):
        self.index+=1
        self.show_question()

    def finish(self):
        try:
            path=quiz_data.save_answers(self.answers,root=self.storage)
        except (OSError,ValueError) as exc:
            messagebox.showerror('Ответы не сохранены',f'{exc}\nОкно останется открытым. Можно повторить сохранение.',parent=self.window)
            return
        if path:
            messagebox.showinfo('Квиз сохранён',f'Ответов: {len(self.answers)}\nФайл для передачи:\n{path}',parent=self.window)
        self.window.destroy()
        self.on_done()
