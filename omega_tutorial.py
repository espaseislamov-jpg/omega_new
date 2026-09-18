"""Operator walkthrough with an isolated practice plot; never edits real data."""
import tkinter as tk
from tkinter import ttk
import os
from pathlib import Path
import numpy as np

from omega_review.quiz_data import outbox

def tutorial_seen_path():
    return Path(os.environ.get('LOCALAPPDATA', str(Path.home()))) / 'Omega' / 'ui' / 'tutorial_v3_seen'


def should_show_tutorial():
    return not tutorial_seen_path().exists()


def mark_tutorial_seen():
    path = tutorial_seen_path()
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.touch(exist_ok=True)
    except OSError:
        import logging
        logging.getLogger('omega').exception('Cannot save tutorial completion')

STEPS = (
    ('Профиль и CSV', 'Начните с профиля прибора',
     'Нажмите «Новый расчёт», выберите профиль своего прибора и откройте CSV. '
     'Программа рассчитает пробы и покажет их слева.\n\n'
     'Через «Настроить профили» можно проверить настройки или создать новый профиль.'),
    ('Проверка серии', 'Начните с проб, требующих внимания',
     'Выбирайте пробу в списке слева или переходите кнопками «Назад» и «Далее». '
     'В фильтре «Требуют обязательной проверки» собраны красные пробы.\n\n'
     'Нажмите на статус справа, чтобы прочитать причины, или на код кислоты под ним — перейти к пику.'),
    ('Чтение графика', 'Смотрите на пик вместе с соседями',
     'Время указано по оси X. Чёрная линия — исходный сигнал после вычитания базовой. '
     'Закрашенная область показывает участок интегрирования.\n\n'
     'Вкладки C16–C22 открывают группы кислот, «Пик» — выбранную кислоту. '
     'Поздние пики можно вернуть флажком «Показать всю запись…».'),
    ('Ручные границы', 'Попробуйте указать две границы',
     'В программе сначала выберите кислоту в таблице, затем нажмите «Выделить границы». '
     'Щёлкните начало и конец пика и нажмите «Применить границы».\n\n'
     'На учебном графике ниже можно потренироваться: выберите две границы и примените их.'),
    ('Добавление пика', 'Новый пик — через правую кнопку мыши',
     'Нажмите ПКМ на графике → «Добавить пик — указать две границы». '
     'Укажите начало и конец, затем выберите кислоту.\n\n'
     'Неназначенный пик виден на графике, но не входит в расчёт Омега-3. '
     'Попробуйте ПКМ на учебном графике ниже.'),
    ('Сохранение', 'Сохраните правки перед передачей работы',
     '«Сохранить правки» выгружает ручные изменения в отдельный файл. '
     '«Загрузить правки» позволяет применить их к соответствующему CSV.\n\n'
     'В «Результатах серии» можно скопировать итоговую таблицу. '
     '«Отмена выделения» снимает выбор кислоты, а «Отменить» отменяет последнюю ручную правку.'),
    ('Квиз и помощь', 'Помогайте улучшать судью',
     'При выходе можно оценить до шести пиков: верно ли назначена кислота, '
     'хороши ли границы, есть ли замечания. Квиз можно пропустить.\n\n'
     'Ответы сохраняются в файл .omega-quiz. Передайте его ответственному за обучение; '
     'он добавит ответы через «Обучение → Импортировать ответы квиза…».'),
)


class OperatorTutorial:
    def __init__(self, app, on_done=None):
        self.app, self.on_done = app, on_done
        self.index, self.closed = 0, False
        self.bounds, self.first, self.inserting = None, None, False
        self.context_popup = None
        self.previous_grab = app.root.grab_current()
        self.window = tk.Toplevel(app.root)
        self.window.title('Как работать с Omega')
        width = min(980, app.root.winfo_screenwidth()-60)
        height = min(720, app.root.winfo_screenheight()-90)
        self.window.geometry(f'{width}x{height}+{max(0,(app.root.winfo_screenwidth()-width)//2)}+{max(0,(app.root.winfo_screenheight()-height)//2)}')
        self.window.minsize(min(820, width), min(590, height))
        self.window.transient(app.root)
        self.window.configure(background='#F3F5F8')
        self.window.protocol('WM_DELETE_WINDOW', self.finish)
        self.window.bind('<Escape>', self.escape)

        head = ttk.Frame(self.window, padding=(22,18,22,12))
        head.pack(fill='x')
        ttk.Label(head, text='Знакомство с Omega', font=('Segoe UI',19,'bold')).pack(anchor='w')
        ttk.Label(head, text='7 коротких шагов · можно пропустить или открыть снова кнопкой «Как работать»',
                  wraplength=750).pack(anchor='w', pady=(5,0))
        bottom = ttk.Frame(self.window, padding=(22,12))
        bottom.pack(side='bottom', fill='x')
        self.back_button = ttk.Button(bottom, text='Назад', command=self.back)
        self.back_button.pack(side='left')
        ttk.Button(bottom, text='Пропустить обучение', command=self.finish).pack(side='left', padx=10)
        self.next_button = ttk.Button(bottom, text='Далее →', style='Accent.TButton', command=self.next)
        self.next_button.pack(side='right')
        self.progress = ttk.Label(bottom)
        self.progress.pack(side='right', padx=18)

        body = ttk.Frame(self.window, padding=(18,0,18,0))
        body.pack(fill='both', expand=True)
        sidebar = ttk.Frame(body)
        sidebar.pack(side='left', fill='y', padx=(0,18))
        self.step_buttons = []
        for index, (label, _, _) in enumerate(STEPS):
            button = ttk.Button(sidebar, text=f'{index+1}. {label}', width=23,
                                command=lambda i=index:self.show(i))
            button.pack(fill='x', pady=4)
            self.step_buttons.append(button)
        self.content = ttk.Frame(body, padding=18, style='Card.TFrame')
        self.content.pack(fill='both', expand=True)
        self.content.bind('<Configure>', self.resize_content)
        self.show(0)
        self.window.grab_set()
        self.next_button.focus_set()

    def resize_content(self, event):
        for widget in self.content.winfo_children():
            if widget.winfo_class() == 'TLabel':
                widget.configure(wraplength=max(220,event.width-40))

    def show(self, index):
        if self.context_popup is not None:
            self.context_popup.destroy()
            self.context_popup = None
        self.index = index
        self.bounds, self.first, self.inserting = None, None, False
        for child in self.content.winfo_children():
            child.destroy()
        for i, button in enumerate(self.step_buttons):
            button.configure(style='Accent.TButton' if i==index else 'TButton')
        self.back_button.state(['disabled'] if index==0 else ['!disabled'])
        self.next_button.configure(text=('Начать работу' if self.on_done else 'Готово') if index==len(STEPS)-1 else 'Далее →')
        self.progress.configure(text=f'{index+1} / {len(STEPS)}')
        _, title, description = STEPS[index]
        ttk.Label(self.content, text=title, font=('Segoe UI',15,'bold'), style='Card.TLabel',
                  wraplength=620).pack(anchor='w', fill='x', pady=(0,12))
        ttk.Label(self.content, text=description, style='Card.TLabel', wraplength=620,
                  justify='left').pack(anchor='w', fill='x')
        self.feedback = tk.StringVar(value='')
        self.feedback_label = ttk.Label(self.content, textvariable=self.feedback, style='CardMuted.TLabel',
                                        wraplength=620, justify='left')
        self.feedback_label.pack(side='bottom', fill='x', pady=(10,0))
        if index == 0:
            self.card('1   Новый расчёт', '2   Профиль вашего прибора', '3   Выбрать CSV и рассчитать')
            self.feedback.set('Новый CSV можно открыть с другим профилем. Эти учебные шаги не меняют настройки программы.')
        elif index == 1:
            for text, bg, fg in [('ПРОВЕРКИ ПРОЙДЕНЫ','#E7F4ED','#25714A'),
                                 ('ПРОВЕРИТЬ — отмеченные пики','#FFF3D8','#805500'),
                                 ('СТОП — обязательная проверка','#FCE4E7','#A52234')]:
                tk.Label(self.content, text=text, background=bg, foreground=fg,
                         anchor='w', font=('Segoe UI',11), padx=14, pady=7).pack(fill='x', pady=(8,0))
            self.feedback.set('Судья помогает выбрать, что проверить. Причины предупреждения доступны по нажатию на статус.')
        elif index in (2,3,4):
            if index in (3,4):
                actions = ttk.Frame(self.content, style='Card.TFrame')
                actions.pack(side='bottom', fill='x', pady=(10,0))
                self.practice_button = ttk.Button(actions, text='Применить границы' if index==3 else 'Добавить пик',
                                                  command=self.apply_practice, style='Accent.TButton')
                self.practice_button.pack(side='left')
                self.practice_button.state(['disabled'])
                ttk.Button(actions, text='Сбросить упражнение', command=lambda:self.show(self.index)).pack(side='right')
                self.feedback.set('Щёлкните две границы на учебном графике.' if index==3 else 'Нажмите правой кнопкой мыши на учебном графике.')
            else:
                self.bounds = (7.574,7.615)
                self.feedback.set('Учебный пример: олеиновая кислота · C18:1N9C. Соседний пик нужен для оценки разделения.')
            self.canvas = tk.Canvas(self.content, background='#FFFFFF', highlightthickness=0, height=220)
            self.canvas.pack(fill='both', expand=True, pady=(16,0))
            self.canvas.bind('<Configure>', lambda event:self.draw_plot())
            self.canvas.bind('<Button-1>', self.click_plot)
            self.canvas.bind('<Button-3>', self.context_menu)
        elif index == 5:
            self.card('CSV + файл правок', 'Результаты серии → копирование таблицы')
            self.feedback.set('После применения границ результат уже пересчитан. Для переноса ручных изменений на другой компьютер сохраните файл правок.')
        else:
            self.card('Как работать → повторить эту обучалку', 'Обучение → Квиз по пикам')
            self.feedback.set(f'Папка файлов квиза на этом компьютере:\n{outbox()}\n\nПока обучалка показывается при каждом запуске.')
        self.window.after_idle(self.refresh_wrap)

    def refresh_wrap(self):
        if not self.closed:
            width = max(220,self.content.winfo_width()-40)
            for widget in self.content.winfo_children():
                if widget.winfo_class() == 'TLabel':
                    widget.configure(wraplength=width)

    def card(self, *lines):
        box = tk.Frame(self.content, background='#EEF4FF', padx=16, pady=8)
        box.pack(fill='x', pady=(12,0))
        for line in lines:
            label = tk.Label(box, text=line, background='#EEF4FF', foreground='#2463C4',
                             font=('Segoe UI',11), anchor='w', justify='left', wraplength=480)
            label.pack(fill='x', pady=3)
            label.bind('<Configure>', lambda event,w=label:w.configure(wraplength=max(180,event.width)))

    def plot_geometry(self):
        return 42., 30., max(100.,self.canvas.winfo_width()-18.), max(55.,self.canvas.winfo_height()-36.)

    def draw_plot(self):
        self.canvas.delete('all')
        left, top, right, bottom = self.plot_geometry()
        x = np.linspace(7.50,7.65,450)
        y = .9*np.exp(-.5*((x-7.555)/.008)**2)+.62*np.exp(-.5*((x-7.592)/.007)**2)
        px = left+(x-7.5)/.15*(right-left)
        py = bottom-y*(bottom-top)*.86
        self.canvas.create_text(left,8,text='УЧЕБНЫЙ ПРИМЕР · C18:1N9C',anchor='nw',fill='#596B80',font=('Segoe UI',9))
        for tick in np.linspace(7.50,7.65,6):
            tx = left+(tick-7.5)/.15*(right-left)
            self.canvas.create_line(tx,top,tx,bottom,fill='#EDF1F6')
            self.canvas.create_text(tx,bottom+13,text=f'{tick:.2f}',fill='#596B80',font=('Segoe UI',8))
        self.canvas.create_text((left+right)/2,bottom+29,text='Время, мин',fill='#596B80',font=('Segoe UI',8))
        self.canvas.create_line(left,bottom,right,bottom,fill='#AAB8CC')
        if self.bounds:
            a,b = self.bounds
            mask = (x>=a)&(x<=b)
            if mask.any():
                points = [px[mask][0],bottom,*np.column_stack((px[mask],py[mask])).ravel(),px[mask][-1],bottom]
                self.canvas.create_polygon(*points,fill='#C9DAF5',outline='')
        self.canvas.create_line(*np.column_stack((px,py)).ravel(),fill='#223047',width=2)
        for time in self.bounds or (() if self.first is None else (self.first,)):
            tx = left+(time-7.5)/.15*(right-left)
            self.canvas.create_line(tx,top,tx,bottom,fill='#2463C4',width=2)

    def context_menu(self, event):
        if self.index != 4:
            return
        if self.context_popup is not None:
            self.context_popup.destroy()
        menu = self.context_popup = tk.Menu(self.window, tearoff=False)
        menu.add_command(label='Добавить пик — указать две границы', command=self.begin_insert)
        try:
            menu.tk_popup(event.x_root,event.y_root)
        finally:
            menu.grab_release()
            if not self.closed:
                self.window.grab_set()

    def begin_insert(self):
        self.inserting = True
        self.first, self.bounds = None, None
        self.practice_button.state(['disabled'])
        self.feedback.set('Теперь щёлкните начало и конец нового пика.')
        self.draw_plot()

    def click_plot(self, event):
        if self.index not in (3,4) or (self.index==4 and not self.inserting):
            return
        left, top, right, bottom = self.plot_geometry()
        if not left<=event.x<=right or not top<=event.y<=bottom:
            return
        time = 7.5+(event.x-left)/(right-left)*.15
        if self.first is None:
            self.first, self.bounds = time, None
            self.practice_button.state(['disabled'])
            self.feedback.set('Первая граница выбрана. Щёлкните вторую.')
        elif abs(time-self.first)<.001:
            self.feedback.set('Границы слишком близко. Щёлкните вторую дальше от первой.')
        else:
            self.bounds = tuple(sorted((self.first,time)))
            self.first = None
            self.practice_button.state(['!disabled'])
            self.feedback.set('Границы готовы. Нажмите кнопку под графиком.')
        self.draw_plot()

    def apply_practice(self):
        if not self.bounds:
            return
        self.practice_button.state(['disabled'])
        self.inserting = False
        if self.index==3:
            if not self.bounds[0]<7.592<self.bounds[1]:
                self.feedback.set('Вы применили границы. Для олеиновой кислоты они должны охватывать правую вершину; попробуйте ещё раз.')
            else:
                self.feedback.set('Готово! В программе после применения пересчитываются площадь и результат. «Отмена выделения» завершает работу с кислотой.')
        else:
            self.feedback.set('Учебный пик добавлен. В программе следующим шагом выберите кислоту или оставьте пик неназначенным.')

    def escape(self, event=None):
        if self.first is not None or self.inserting:
            self.show(self.index)
        else:
            self.finish()
        return 'break'

    def back(self):
        if self.index:
            self.show(self.index-1)

    def next(self):
        if self.index==len(STEPS)-1:
            self.finish()
        else:
            self.show(self.index+1)

    def finish(self):
        if self.closed:
            return
        self.closed = True
        self.window.grab_release()
        self.window.destroy()
        self.app.tutorial_window = None
        if self.previous_grab is not None and self.previous_grab.winfo_exists():
            self.previous_grab.grab_set()
        if self.on_done is not None and self.app.root.winfo_exists():
            self.on_done()
