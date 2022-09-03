
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
#-----------------------------------------------Исходные данные--------------------------------------------------------
time_window_us = 300     # временной интервал на котором рассматривается корреляция, мкс
sigma_normal = 11          # СКО гауссовского случайного процесса
mu_normal = 0
corr_int_us = 11         # интеравал корреляции, мкс
n_samples = 20              # число усредняемых реализаций
f_sampling_mhz = 20         # частота дискретизации, МГц
n_counts = np.ceil(time_window_us*f_sampling_mhz); # число отсчетов
#----------------------------------------------------------------------------------------------------------------------
edges = 37


def corr_filter(Z, rL_us, corr_us):
    signal_length = Z.shape[1]
    x = np.linspace(-rL_us/2,rL_us/2,signal_length)
    F = np.exp(-abs(x) / (corr_us / 2))
    f = np.real(np.sqrt(2)*np.fft.ifft(np.fft.fft(Z,signal_length) * np.sqrt(np.fft.fft(F,signal_length)),signal_length))
    return f


def corr_filter_linear(Z, rL_us, corr_us):
    signal_length = Z.shape[1]

    x = np.linspace(-rL_us / 2, rL_us / 2, signal_length)
    F = np.ndarray(shape=(signal_length))

    for i in range(signal_length):
        if abs(x[i]) < corr_us:
            F[i] = 1 - abs(x[i]) / corr_us
        else:
            F[i] = 0

    f = np.real(np.sqrt(2) * np.fft.ifft(np.fft.fft(Z, signal_length) * np.sqrt(np.fft.fft(F, signal_length)), signal_length))
    return f


def corrf_lin(x,corr_int_us):
    anal_cor_fun = np.ndarray(shape=len(x))
    for i in range(len(x)):
        if abs(x[i]) < corr_int_us:
            anal_cor_fun[i] = 1 - abs(x[i]) / corr_int_us
        else:
            anal_cor_fun[i] = 0
    return  anal_cor_fun


def one_realiz(t,y):
    r_out_single = np.correlate(y, y, 'full')
    shift_axis_us = np.linspace(-int(time_window_us), int(time_window_us), 2 * int(n_counts) - 1)

    #---plotting------------------
    fig, (ax1, ax2) = plt.subplots(2,1)
    ax1.plot(t, y, 'b-')
    ax1.set_title('Исходная реализация')
    #ax1.legend(loc='best', frameon=False)
    ax1.set_xlabel(r'$\mu$'+'s')
    ax1.set_ylabel('A',rotation='horizontal')

    ax2.plot(shift_axis_us,r_out_single)
    ax2.set_title('АКФ исходной реализации',fontsize=10)
    ax2.set_xlabel('f KHz')
    ax2.set_ylabel('G(f)',rotation='horizontal',loc='top')
    plt.tight_layout()
    plt.show()


def plot_corr (shift_axis_us, rout, rout_anal):

    plt.figure()
    plt.plot(shift_axis_us,rout)
    plt.plot(shift_axis_us, rout_anal)
    plt.show()


def gistogram_plot(X,analytic,t, title):
    fig, ax = plt.subplots(1, 1, num = 'Нормализированная гистограмма исходного процесса')

    ax.plot(t, analytic, 'g-', label='analytic')

    ax.hist(X, density=True, bins=t, weights = np.ones_like(X)/len(X), label = 'КФ- ' + title )

    ax.legend(loc='best', frameon=False)
    plt.xlabel('x')
    plt.ylabel('W(x)',rotation='horizontal', loc='top')
    plt.title('Нормализированная гистограмма исходного процесса')
    plt.show()


def gistogram_2plot(X1,X2,t, title):
    fig, ax = plt.subplots(1, 1, num='Нормализированная гистограмма КФ-' + title)

    ax.hist([X1, X2], density=True, bins=t,  label = ['Исходный процесс' ,'КФ- ' + title] ) #weights = np.ones_like(X)/len(X),

    ax.legend(loc='best', frameon=False)
    plt.xlabel('x')
    plt.ylabel('W(x)',rotation='horizontal', loc='top')
    plt.title('Нормализированная гистограмма КФ-' + title)
    plt.show()


def mean_realization(y):
    rout = np.zeros((int(n_samples), int(2 * n_counts - 1)))  # задание массива array
    for i in range(1, n_samples - 1, 1):
        rout[i, :] = np.correlate(y[i], y[i], 'full')
    rout_mean = np.mean(rout, axis=0)  # усредненная реализация
    return rout_mean


def plot_SPM(shift_axis_us,X,title):

    spm = np.zeros((int(n_samples), int(2 * n_counts - 1)))  # задание массива array
    for i in range(1, n_samples - 1, 1):
        spm[i, :] = np.fft.fftshift(abs(np.fft.fft(np.correlate(X[i], X[i], 'full'))))
    spm_mean = np.mean(spm, axis=0)  # усредненная реализация

    F = np.linspace(10 ** 3/min(shift_axis_us), 10 ** 3/max(shift_axis_us), shift_axis_us.size)

    plt.figure(num='SPM ' + title)
    plt.plot(F , spm_mean)#, label = 'СПМ ' + title)
    plt.xlabel('f KHz')
    plt.ylabel('G(f)',rotation='horizontal', loc='top')
    plt.title('СПМ ' + title)
    plt.show()


 # ---------------------построение на одном графике -----------------------
def plot_KF(y,y2,anal_cor_fun1,anal_cor_fun2,shift_axis_us,title):
    plt.figure(num='KF '+title)

    plt.plot(shift_axis_us, mean_realization(y) / max(mean_realization(y)), label='y_mean 1 cor int ')
    plt.plot(shift_axis_us, mean_realization(y2) / max(mean_realization(y2)), label='y_mean 2 cor int ')
    plt.plot(shift_axis_us, anal_cor_fun1, label='analytic 1 cor int')
    plt.plot(shift_axis_us, anal_cor_fun2, label='analytic 2 cor int')
    plt.legend(loc='best', frameon=False, title = 'KF- ' + title)

    plt.xlabel(r'$\tau$')
    plt.ylabel(r'R($\tau$)',rotation='horizontal', loc='top')
    plt.title('Усредненные корреляционные функции КФ- ' + title)
    plt.show()


shift_axis_us = np.linspace(-int(time_window_us), int(time_window_us), 2 * int(n_counts) - 1)
#------------------построение исходной реализации-----------------------
x = np.random.normal( loc=mu_normal, scale=sigma_normal, size=[int(n_samples),int(n_counts)])

AXIS = np.linspace(-3*sigma_normal + mu_normal, 3*sigma_normal + mu_normal,edges)
centers = 0.5 * (AXIS[1:]+AXIS[:-1])

gistogram_plot(x[1],norm.pdf(centers,loc=mu_normal, scale=sigma_normal), centers, title= 'Отсутсвует')

t = np.linspace(0,int(time_window_us),int(n_counts))
one_realiz(t,x[0])
plot_SPM(shift_axis_us, x, 'Исходной реализации')

#------------------построение реализации W(x')---------------------------


y = corr_filter(x, time_window_us, corr_int_us)
y2 = corr_filter(x, time_window_us, 2 * corr_int_us) #удвоенное время кор

anal_cor_fun1 = np.exp(-abs(shift_axis_us) / (corr_int_us / 2))
anal_cor_fun2 = np.exp(-abs(shift_axis_us) / (2 * corr_int_us / 2))

gistogram_2plot(y[1],x[1], centers, title= 'экспоненциальная' ) #norm.pdf(centers,loc=mu_normal, scale=sigma_normal)

#one_realiz(t,y[0])
plot_SPM(shift_axis_us, y, 'КФ - экспоненциальная')
plot_KF(y,y2,anal_cor_fun1,anal_cor_fun2,shift_axis_us,title= 'экспоненциальная')

#------------------построение реализации W(x')но с другой КФ---------------------------
y3 = corr_filter_linear(x, time_window_us, corr_int_us)
y4 = corr_filter_linear(x, time_window_us, 2 * corr_int_us) #удвоенное время кор

gistogram_2plot(y3[1],x[1], centers, title= 'линейная' )

anal_cor_fun3 = corrf_lin(shift_axis_us,corr_int_us)
anal_cor_fun4 = corrf_lin(shift_axis_us,corr_int_us * 2)

plot_KF(y3,y4,anal_cor_fun3,anal_cor_fun4,shift_axis_us,title= 'линейная')
plot_SPM(shift_axis_us, y3, 'КФ - линейная')