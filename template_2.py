#%% Смоделировать гауссовский случайный процесс с заданной КФ
"""
В программе: 
1) Загрузка дополнительных библиотек
2) Определение функции фильтрации, обеспечивающей заданную КФ
3) Инициализация параметров модели
4) Формирование некоррелированного СП x
5) Получение из него коррелированного СП y с заданной функцией корреляции
6) Оценка КФ полученного СП y по его реализациям
7) Построение графиков, содержащих 
    - реализации процессов
    - реализации КФ
    - оценочная КФ, полученная усреднением множества КФ и аналитически заданная КФ.

Для вывода графиков в отдельное окно используйте
>>> %matplotlib qt
"""
#%% Импорт библиотек и определения функций
""" Импорт библиотек и определения функций
"""
import numpy as np
from numpy.fft import rfft, irfft

import matplotlib.pyplot as plt
sqrt = np.sqrt


# Фильтрация через fft
def cyclic_concvolution(sig, ir):
    """  
    Циклическая свёртка сигнала с импульсной характеристикой (ИХ) через спектральный метод.

    Функция выполняет свёртку во временной области с использованием преобразования Фурье.
    Для корректной работы требуется равенство длин сигнала и импульсной характеристики.

    Алгоритм:
    1. Проверяется равенство длин входного сигнала и ИХ.
    2. Вычисляется частотная характеристика фильтра как квадратный корень из БПФ ИХ.
    3. Частотная характеристика нормируется на максимальное значение, чтобы избежать
       усиления сигнала.
    4. Выполняется обратное преобразование Фурье после умножения спектров.
    5. Результат возвращается в действительной форме.

    Параметры
    ----------
    sig : np.ndarray
        Входной сигнал (массив чисел).
    ir : np.ndarray
        Импульсная характеристика фильтра, длина которой должна совпадать с длиной сигнала.

    Возвращает
    -------
    np.ndarray
        Отфильтрованный сигнал после циклической свёртки.

    Исключения
    ----------
    ValueError
        Если длина сигнала и импульсной характеристики не совпадают.

    Пример
    -------
    >>> import numpy as np
    >>> from numpy.fft import rfft, irfft
    >>> sig = np.random.randn(128)
    >>> ir = np.hanning(128)
    >>> y = cyclic_concvolution(sig, ir)
    """

    if sig.shape[-1]==ir.shape[-1]:     length = sig.shape[-1]
    else:   raise ValueError('Для этого метода длина ИХ должна быть равна длине сигнала')
    freq_resp = sqrt(rfft(ir,length))    # АЧХ фильтра
    freq_resp = freq_resp/np.max(np.abs(freq_resp))                                 # нормировка АЧХ на максимальное значение: частоты в не должны усиливаться
    sig = np.real(np.sqrt(2)*irfft(rfft(sig, length) * freq_resp,length))
    return sig


# Вывод результатов
def panel_realization(in_signal, out_signal, Rx, Ry, time_window_us, cor_lims=None):
    """
    Построение панели графиков для анализа одной реализации случайного процесса
    и его коррелированной версии.

    Функция формирует фигуру с двумя подграфиками:
    1. Сравнение одной реализации исходного сигнала и сигнала после корреляционной обработки.
    2. Корреляционные функции (КФ) этих сигналов для той же реализации.

    Алгоритм
    --------
    - Строится временная ось реализации (`t_axis`) в пределах указанного окна.
    - Строится ось сдвига (`shift_axis_us`) для отображения корреляционной функции.
    - В первом подграфике отображаются временные зависимости входного и выходного сигнала.
    - Во втором подграфике отображаются корреляционные функции исходного и коррелированного сигналов.
    - Для обоих графиков добавляются сетка, подписи осей и легенда.

    Параметры
    ----------
    in_signal : np.ndarray
        Двумерный массив входного сигнала, где ось 0 — это номер реализации, ось 1 — отсчёты.
    out_signal : np.ndarray
        Двумерный массив выходного (коррелированного) сигнала той же размерности, что и `in_signal`.
    Rx : np.ndarray
        Корреляционная функция исходного сигнала (одной реализации).
    Ry : np.ndarray
        Корреляционная функция выходного (коррелированного) сигнала (той же реализации).
    time_window_us : float
        Длительность временного окна в микросекундах, по которому строится реализация сигнала.

    Возвращает
    ----------
    None
        Функция выполняет построение графиков и не возвращает значения.

    Пример
    ----------
    >>> import numpy as np
    >>> in_signal = np.random.randn(10, 256)
    >>> out_signal = np.random.randn(10, 256)
    >>> Rx = np.random.randn(2*256-1)
    >>> Ry = np.random.randn(2*256-1)
    >>> panel_realization(in_signal, out_signal, Rx, Ry, time_window_us=50)
    """
    n_counts = in_signal.shape[-1]
    t_axis = np.linspace(0, time_window_us, n_counts)                           # ось времени реализации сигнала
    shift_axis_us = np.linspace(-time_window_us, time_window_us, 2*n_counts-1)  # ось отстройки КФ
    if cor_lims is None: cor_lims = [-time_window_us/2, time_window_us/2]
    
    plt.figure()
    plt.subplot(121)                                    # Одна реализация
    plt.plot(t_axis, in_signal[1,:], t_axis, out_signal[1,:])
    plt.title('Одна реализация СП')
    plt.legend(('Исходный','Коррелированный'))
    plt.xlabel('t'), plt.ylabel('U(t)', rotation='horizontal')
    plt.grid(True)

    plt.subplot(122)                                    # КФ одной реализации
    plt.plot(shift_axis_us, Rx, label='Исходный',)
    plt.plot(shift_axis_us, Ry, label='Коррелированный')
    plt.xlim(-3*corr_int_us, 3*corr_int_us)
    plt.title('КФ одной реализации')
    plt.legend()
    plt.xlabel(r'$\tau$'), plt.ylabel(r'R($\tau$)', rotation='horizontal')
    plt.grid(True)



#%% Параметры модели 
""" Параметры модели 
"""
time_window_us = 300    # временной интервал на котором рассматривается корреляция, мкс
f_sampling_mhz = 20     # частота дискретизации, МГц
sigma_normal = 11       # СКО гауссовского случайного процесса
mu_param = 0            # мат ожидание исходного процесса
corr_int_us = 10        # интеравал корреляции, мкс
n_samples = 10          # число усредняемых реализаций
n_counts = int(time_window_us*f_sampling_mhz); # число отсчетов


#%% Определение АКФ 
""" Определение АКФ 
"""
tau_ax = np.arange(-time_window_us/2,time_window_us/2,1/f_sampling_mhz)     # шкала времени для АКФ
target_acf = np.exp(-abs(tau_ax) / (corr_int_us / 2))                       # целевая корреляционная функция

plt.figure()
plt.plot(tau_ax, target_acf, label='Целевая')
plt.legend()
plt.title('Целевая АКФ')
plt.xlabel('Δτ, мкс')
plt.grid(True)


#%% Моделирование случайных процессов
""" Моделирование случайных процессов
"""
x = np.random.normal(mu_param, sigma_normal, (n_samples,n_counts))  # исходный СП
y = cyclic_concvolution(x, target_acf)                              # коррелированный СП


#%% Оценка КФ
""" Оценка КФ
- Рассчитать
"""
rx_sample = np.correlate(x[1,:], x[1,:], mode='full') / np.var(x[1,:]) / n_counts
ry_est = np.zeros((n_samples, 2*n_counts-1))
for i in range(n_samples):          # построчный расчёт корреляции
    ry_est[i,:]=np.correlate(y[i,:], y[i,:], mode='full')/ np.var(y[i,:]) / n_counts     # КФ
ry_sample = ry_est[1,:]                                                                  # КФ одной шумовой рализации
ry_est = np.mean(ry_est, axis=0)                                                          # усреднение КФ по множеству реализаций


#%% Вывод результатов 
""" Вывод результатов
"""
panel_realization(x,y, rx_sample, ry_sample, time_window_us)

shift_axis_us = np.linspace(-time_window_us, time_window_us, 2*n_counts-1)  # ось отстройки КФ
r_fun = np.exp(-abs(shift_axis_us) / (corr_int_us / 2))     # аналитически заданная КФ

plt.figure()
plt.subplot(2,2,(3,4))                              # Усреднённая КФ
plt.plot(shift_axis_us, ry_est, shift_axis_us, r_fun)
plt.xlim(-3*corr_int_us, 3*corr_int_us)
plt.title('Усреднённая КФ')
plt.legend(('моделирование','аналитически заданная'))
plt.xlabel(r'$\tau$'), plt.ylabel(r'R($\tau$)', rotation='horizontal');
plt.grid(True)
plt.show()


# %%
# Выполните моделирование СП с КФ заданной по варианту
# Проанализируйте распределение отсчётов СП до и после фильтра

plt.figure()
plt.hist([x.flatten(),y.flatten()], bins=64, density=True, label = ['исх','результ'])
plt.legend(loc='best', frameon=False)
plt.title('Нормальное распределение, m=' + str(mu_param) + ', sigma= ' + str(sigma_normal))
plt.xlabel('y')
plt.ylabel('W(y)', rotation='horizontal')
plt.show()