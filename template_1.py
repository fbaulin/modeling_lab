#%% Сформировать нормально распроеделенную СВ на основе равномерно распределенных СВ
"""
В программе:
1) Загрузка дополнительных библиотек
2) Инициализация параметров моделирования и распределений
3) Формирование СВ с использованем встроенной функции
4) Формировнаие СВ путём комбинирования Релеевской и арксинусной СВ
5) Задание аналитической функции распределения
6) Построение графика, содержащего:
    - аналитически заданную ПВ
    - оценённую по СВ, полученной с использованием внутренней функции
    - оценённую по СВ, полученной на основании преобразований.

Для вывода графиков в отдельное окно используйте
>>> %matplotlib qt
"""
#%% Загрузка дополнительных библиотек
""" Загрузка дополнительных библиотек и определние 
"""
# from scipy.stats import norm, uniform, expon
from numpy.random import normal, uniform
import matplotlib.pyplot as plt
import numpy as np

# Вывод панели [сигнал | распределение]
def signal_pdf_panel(labels, signals, bins=None, ideal_pdf_w=None, ideal_pdf_x=None, limit=200, fig=None):
    """ 
    Построение графиков анализа случайной величины.

    Функция отображает два графика:
    1. Реализацию сигнала (первые `limit` отсчётов временного ряда).
    2. Гистограмму распределения выборки (или нескольких выборок) с 
    наложенной идеализированной плотностью распределения вероятностей (ПРВ).

    Параметры
    ----------
    labels : list of str
        Подписи для легенды (названия методов генерации).
    signals : list of array_like 
        Список последовательностей случайных величин (выборки).
    bins : int, optional
        Количество интервалов в гистограмме. По умолчанию 10.
    ideal_pdf_w : array_like, optional
        Значения случайной величины (ось X для правого графика). 
        Если None, то идеализированная плотность не отображается.
    ideal_pdf_x : array_like, optional
        Значения плотности вероятности, соответствующие `ideal_pdf_w`. 
        Если None, то идеализированная плотность не отображается.
    limit : int, optional
        Длина фрагмента реализации, который строится в левом окне. По умолчанию 200.
    fig : matplotlib.figure.Figure, optional
        Существующая фигура для добавления графиков. 
        Если None — создаётся новая.

    Возвращает
    -------
    fig : matplotlib.figure.Figure
        Фигура Matplotlib с построенными графиками.

    Пример
    -------
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> from scipy.stats import uniform
    >>> data1 = np.random.rand(1000)
    >>> data2 = uniform.rvs(size=1000)
    >>> fig = signal_pdf_panel(
    ...     labels=["numpy.rand", "scipy.uniform"], 
    ...     signals=[data1, data2], 
    ...     bins=20, 
    ...     limit=100
    ... )
    
    """
    if fig is None:     # создать новое оккно
        fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    else:               # использовать существующее окно в fig
        axes = fig.get_axes()

    # Left: Time series comparison (first 100 points)
    if isinstance(signals,list):
        for s, l in zip(signals,labels):
            axes[0].plot(s[:limit], label=l, lw=1)
    else:   axes[0].plot(signals[:limit], label=labels, lw=1)
    axes[0].set_title("Реализация сигнала")
    axes[0].set_xlabel("Время")
    axes[0].set_ylabel("U, В")
    axes[0].legend()
    axes[0].grid(True, alpha=0.1, color='black', zorder=0)

    # Right: Histograms + ideal PDF
    if ideal_pdf_w is not None and ideal_pdf_x is not None:
        axes[1].plot(ideal_pdf_w, ideal_pdf_x, 'r-', lw=1.5, label="Аналитическая")
    # axes[1].hist(signals, bins=bins, density=True, orientation='horizontal',
    #             alpha=0.5, label=labels, edgecolor='Grey')
    axes[1].hist(signals, bins=bins, density=True, orientation='horizontal',
                label=labels, edgecolor=None)
    axes[1].set_title("Плотности РВ")
    axes[1].set_xlabel("W(U)")
    axes[1].legend()
    axes[1].grid(True, alpha=0.1, color='black', zorder=0)

    plt.tight_layout()
    return fig


# Вывод панели распредедений
def pdf_panel(labels, signals, bins=None, ideal_pdf_w=None, ideal_pdf_x=None, fig=None, title=None):
    """
    Построение гистограмм распределения случайных величин и, при наличии, 
    наложение аналитической плотности вероятности.

    Функция отображает одну ось с гистограммами (одной или нескольких выборок) 
    и при передаче аргументов `ideal_pdf_w` и `ideal_pdf_x` строит 
    идеализированную (аналитическую) кривую плотности.

    Параметры
    ----------
    labels : list of str
        Список подписей для легенды, соответствующих выборкам из `signals`.
    signals : list of array_like
        Список выборок случайных величин, для которых строятся гистограммы.
    bins : int, optional
        Количество интервалов в гистограмме. По умолчанию используется 
        стандартное значение matplotlib.
    ideal_pdf_w : array_like, optional
        Значения плотности вероятности (ось Y).
    ideal_pdf_x : array_like, optional
        Значения случайной величины, соответствующие `ideal_pdf_w` (ось X).
    fig : matplotlib.figure.Figure, optional
        Существующая фигура для добавления графиков. Если None, создаётся новая.
    title : str, optional
        Заголовок графика.

    Возвращает
    -------
    ax : matplotlib.axes.Axes
        Ось Matplotlib с построенными графиками.

    Пример
    -------
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>>
    >>> data1 = np.random.rand(1000)
    >>> data2 = np.random.rand(1000)
    >>> x = np.linspace(0, 1, 200)
    >>> pdf = np.ones_like(x)  # аналитическая плотность для U[0,1]
    >>>
    >>> ax = pdf_panel(
    ...     labels=["Генератор 1", "Генератор 2"],
    ...     signals=[data1, data2],
    ...     bins=20,
    ...     ideal_pdf_w=pdf,
    ...     ideal_pdf_x=x,
    ...     title="Сравнение распределений"
    ... )
    >>> plt.show()
    """


    if fig is None:     plt.figure()
    else:               plt.sca(fig.get_axes())
    ax = plt.gca()
    if (ideal_pdf_w is not None) and (ideal_pdf_x is not None):
        plt.plot(ideal_pdf_x,ideal_pdf_w,'r-',label='Аналитичекая')
    ax.hist(signals, density=True, bins=bins, label=labels)
    plt.legend(loc='best', frameon=False)
    plt.title(title)
    plt.xlabel('y')
    plt.ylabel('W(y)', rotation='horizontal')
    return fig


#%% Моделирование равномерной СВ
""" Моделирование равномерной СВ
Задачи:
1.  Сгенерировать равномерную СВ, используя пакет scipy, и построить 
    последовательность СВ и распределение (вид гистограммы и идеализированный вид).
2.  Реализовать конгруэнтный генератор равномерной СВ, распределенной 
    в интервале [0; 1) и построить последовательность и плотность распределения.
3.  Построить на одном графике гистограммы полученные в результате использования 
    библиотечного метода и своего генератора,
    а также идеализированную плотность распределения.

Библиотечная функция, используемая для формирования равномерной СВ 
np.random.uniform или scipy.stats.uniform.
Примеры использования функции построения графиков
Один процесс
>>> signal_pdf_panel('Процесс X', x ,ideal_pdf_w=uniform_ideal_w, ideal_pdf_x=uniform_ideal_x)
Два процесса
>>> signal_pdf_panel(['Процесс X', 'Процесс Y'], [x, y] ,ideal_pdf_w=uniform_ideal_w, ideal_pdf_x=uniform_ideal_x)
"""
n_values=1000

uniform_ideal_x = np.arange(-0.5,1.5,0.01)
uniform_ideal_w = np.where((uniform_ideal_x>=0)&(uniform_ideal_x<1),1.0,0)

x = uniform(0,1,n_values)

signal_pdf_panel('Библиотечная', x ,
                 ideal_pdf_w=uniform_ideal_w, ideal_pdf_x=uniform_ideal_x);

# Конгруэнтный метод


#%% Моделирование нормальной СВ
""" Моделирование нормальной СВ
Задачи:
1.  Сгенерировать выборку БГШ, используя пакет scipy, и построить
    фрагмент этой выборки и распределение (вид гистограммы и идеализированный вид).
2.  Сгенерировать БГШ методом произведения релеевской и арксинусной СВ (см. конспект), 
    на основании выборки равномерных СВ,сгенерированных библиотечной функцией 
    и построить выборку полученного БГШ и плотность распределения.
3.  Построить на одном графике гистограммы полученные в результате использования 
    библиотечного метода и своего генератора,
    а также идеализированную плотность распределения.
4.  Вывести плотности распределения исходного процесса и результирующего процесса.

Для построения графиков можете использовать функцию из предыдущей секции.
Для формирования нормальной величины используйте np.random.normal.
"""
# Параметры моделирования и распределения
sigma_param = 1
mu_param=0          
n_values=1000
# Формирование нормально распределенных СВ
y_lib = normal(loc=mu_param,scale=sigma_param,size=n_values) # встроенная функция нормального распределения

# Формиронование нормально распределенной СВ (см. раздел 1.2.7 в конспекте)
x_1 = uniform(0,1,size=n_values)      # исходные равномерные СВ  
x_2 = uniform(0,1,size=n_values)      # исходные равномерные СВ

y_1 = (-2 * np.log(x_1)) ** (1/2)       # результат н.лин. преобразования - Релеевская СВ
y_2 = np.sin(2*np.pi*x_2)               # результат н.лин. преобразования - Арксинусная СВ
y_custom = sigma_param * y_1 * y_2 + mu_param   # комбинирование для получения нормальной СВ

#%% Вывод ПВ исходной равномерной СВ
n_edges = 20
x_hist = np.linspace(-3*sigma_param + mu_param, 3*sigma_param + mu_param,n_edges)
uniform_ideal_x = np.arange(-0.5,1.5,0.01)
uniform_ideal_w = np.where((uniform_ideal_x>=0)&(uniform_ideal_x<1),1.0,0)

pdf_panel(['X1', 'X2'], [x_1, x_2], 
          ideal_pdf_w=uniform_ideal_w, ideal_pdf_x=uniform_ideal_x, 
          title='Исходное распределение');

#%% Вывод ПВ результирующей нормальной СВ


# %% 
# Моделирование распределения заданного по варианту 
# Скопируйте текст программы в первом разделе и модифицируйте его так, чтобы получить нужное распределение

