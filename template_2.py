# %% Смоделировать гауссовский случайный процесс с заданной КФ
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
# Загрузка дополнительных библиотек
import numpy as np
import matplotlib.pyplot as plt


# Фильтрация для получения образцовой КФ e^(-2*abs(x)/corr)
def corr_filter(Z, rL_us, corr_us):
    signal_length = Z.shape[1]
    x = np.linspace(-rL_us/2,rL_us/2,signal_length)
    f = np.exp(-abs(x) / (corr_us / 2))             # целевая корреляционная функция
    F = np.sqrt(np.fft.fft(f,signal_length))        # АЧХ фильтра
    F = F/np.max(np.abs(F))                         # нормировка АЧХ на максимальное значение: гармоники не должны училиваться
    f = np.real(np.sqrt(2)*np.fft.ifft(np.fft.fft(Z,signal_length) * F,signal_length))
    return f


# %%
# Формирование реализаций и расчёт
# Параметры моделирования и распределения
time_window_us = 300    # временной интервал на котором рассматривается корреляция, мкс
f_sampling_mhz = 20     # частота дискретизации, МГц
sigma_normal = 11       # СКО гауссовского случайного процесса
mu_param = 0            # мат ожидание исходного процесса
corr_int_us = 10        # интеравал корреляции, мкс
n_samples = 20          # число усредняемых реализаций
n_counts = int(time_window_us*f_sampling_mhz); # число отсчетов

# Моделирование случайных процессов
x = np.random.normal(mu_param, sigma_normal, (n_samples,n_counts))  # исходный СП
y = corr_filter(x, time_window_us, corr_int_us)                     # коррелированный СП
# Расчёт КФ
r_est = np.zeros((n_samples, 2*n_counts-1))
for i in range(n_samples):          # построчный расчёт корреляции
    r_est[i,:]=np.correlate(y[i,:], y[i,:], mode='full')/ np.var(y[i,:]) / n_counts    # КФ
ry_sample = r_est[1,:]               # КФ одной шумовой рализации
r_est = np.mean(r_est, axis=0)      # усреднение КФ по множеству реализаций
rx_sample = np.correlate(x[1,:], x[1,:], mode='full') / np.var(x[1,:]) / n_counts

# %% 
# Графика
t_axis = np.linspace(0, time_window_us, n_counts)                           # ось времени реализации сигнала
shift_axis_us = np.linspace(-time_window_us, time_window_us, 2*n_counts-1)  # ось отстройки КФ

r_fun = np.exp(-abs(shift_axis_us) / (corr_int_us / 2))     # аналитически заданная КФ

plt.figure()
plt.subplot(221)                                    # Одна реализация
plt.plot(t_axis, x[1,:], t_axis, y[1,:])
plt.xlim(0, 10*corr_int_us)
plt.title('Одна реализация СП')
plt.legend(('исходный','коррелированный процесс'))
plt.xlabel('t'), plt.ylabel('U(t)', rotation='horizontal')
plt.subplot(222)                                    # КФ одной реализации
plt.plot(shift_axis_us, rx_sample, shift_axis_us, ry_sample)
plt.xlim(-3*corr_int_us, 3*corr_int_us)
plt.title('КФ одной реализации')
plt.legend(('исходный','коррелированный процесс'))
plt.xlabel(r'$\tau$'), plt.ylabel(r'R($\tau$)', rotation='horizontal')
plt.subplot(2,2,(3,4))                              # Усреднённая КФ
plt.plot(shift_axis_us, r_est, shift_axis_us, r_fun)
plt.xlim(-3*corr_int_us, 3*corr_int_us)
plt.title('КФ')
plt.legend(('моделирование','аналитически заданная'))
plt.xlabel(r'$\tau$'), plt.ylabel(r'R($\tau$)', rotation='horizontal');
plt.show()


# %%
# Выполните моделирование СП с КФ заданной по варианту
# Проанализируйте распределение отсчётов СП до и после фильтра

plt.figure()
plt.hist([x.flatten(),y.flatten()], bins=64, density=True, label = ['исх','результ'])
plt.legend(loc='best', frameon=False)
plt.title('Нормальное распределение, m=' + str(mu_param) + ', sigma= ' + str(sigma_param))
plt.xlabel('y')
plt.ylabel('W(y)', rotation='horizontal')
plt.show()