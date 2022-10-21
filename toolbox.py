import matplotlib.pyplot as plt
import math
import numpy as np
from scipy.fft import fft, fftfreq, fftshift, ifft, rfft
from scipy.signal import butter, lfilter, cheby2, cheby1, impulse2, freqz, spectrogram
import pywt


# Отображение сигналов
def plot_signal(signal_args):
#PLOT_SIGNAL Отобразить сигнал, в случае нескольких сигналов, они
#передаются по горизонтали
# Аргументы:
#   один аргумент - значения сигнала, тогда по оси абсцис отсчеты
#   два аргумента - {временная шкала или шаг дискретизации, значения амплитуд}
#   три аргумента - третьим передается название окна
    lines = ["-","--","-.",":"]
    leg = []
    plt.figure()
    for i in range (len(signal_args)):
        markers = {'-r','--g',':b'}
        args_len = len(signal_args[0])
        if (args_len > 3 or args_len < 1):
            print('Неверный формат ввода')
            return
        if (args_len == 1):
            signal = signal_args[i][0]
            x_axis = range(0,len(signal))
        if (args_len == 2 or args_len == 3):
            signal = signal_args[i][1]
            if (type(signal_args[i][0]) == float):
                x_axis = np.arange(0,signal_args[i][0]*(len(signal)),signal_args[i][0])
            else:
                x_axis = signal_args[i][1]
            if (args_len == 2):
                tab_name ='Сигнал'
            else:
                tab_name = signal_args[i][2]
                
    
         # Построение графика
        
        plt.title("Name") # заголовок
        plt.xlabel("t, мкc") # ось абсцисс
        
        plt.plot(x_axis*1e6, signal, linestyle=lines[i])  # построение графика
        
        if (args_len == 3):
            leg.append(signal_args[i][2])
        
    plt.legend(leg)
    plt.title('Исходный сигнал')
    plt.show()


# Отображение спектр сигнала
def plot_spectum(signal_array):
    t_d_us = signal_array[0][0]*1e6
    #plt.figure()
    fig, axs = plt.subplots(len(signal_array), 1)
    
    for i in range(len(signal_array)):
        signal = signal_array[i][1]
        n_signal = signal.size
        t_window = n_signal*t_d_us
        f_step_mhz = 1/t_window
        if n_signal % 2 == 0:
            f_axis = np.arange(0,(n_signal+1)*f_step_mhz/2,f_step_mhz)
            tmp_axis = np.flip(f_axis[2:]) * (-1)
            f_axis = np.append(tmp_axis, f_axis)
        else:
            print('нечетное число отсчетов')
            f_axis = np.arange(0,(n_signal)*f_step_mhz/2,f_step_mhz)
            tmp_axis = np.flip(f_axis[1:]) * (-1)
            f_axis = np.append(tmp_axis, f_axis)
    
        signal_spectrum = 1/n_signal*fft(signal, n_signal)
        if len(signal_array) == 1:
            plt.title("Спектр") # заголовок
            plt.xlabel("Частота, МГц") # ось абсцисс
        
            plt.plot(f_axis, fftshift(abs(signal_spectrum)))  # построение графика
        else:
            fig.suptitle('Спектр')
             # Построение графика
            if (len(signal_array[i]) > 2):
                axs[i].set(label = (signal_array[i][2]))
        
            # # заголовок
            axs[i].set(xlabel = ("Частота, МГц")) # ось абсцисс
            axs[i].plot(f_axis, fftshift(abs(signal_spectrum)), label=(signal_array[i][2]))  # построение графика
            axs[i].legend()
    plt.show()


# Отображение спектрограммы
def plot_spectrograms(data):
    plt.figure()
    plt.subplot(311)

    #большое количество точек для БПФ
    powerSpectrum, freqenciesFound, time, imageAxis = plt.specgram(data, NFFT=1024, Fs=None, Fc=None, detrend=None, window=None, 
             noverlap=None, cmap=None, xextent=None, pad_to=None, sides=None, 
             scale_by_freq=None, mode=None, scale=None, vmin=None, vmax=None)

    plt.xlabel('Sample')
    plt.ylabel('Normalized Frequency')
    plt.subplot(312)

    #стандартный шаг
    powerSpectrum, freqenciesFound, time, imageAxis = plt.specgram(data, NFFT=None, Fs=None, Fc=None, detrend=None, window=None, 
             noverlap=None, cmap=None, xextent=None, pad_to=None, sides=None, 
             scale_by_freq=None, mode=None, scale=None, vmin=None, vmax=None)

    plt.xlabel('Sample')
    plt.ylabel('Normalized Frequency')
    plt.subplot(313)

    #маленький шаг
    powerSpectrum, freqenciesFound, time, imageAxis = plt.specgram(data, NFFT=None, Fs=None, Fc=None, detrend=None, window=None, 
             noverlap=0, cmap=None, xextent=None, pad_to=None, sides=None, 
             scale_by_freq=None, mode=None, scale=None, vmin=None, vmax=None)

    plt.xlabel('Sample')
    plt.ylabel('Normalized Frequency')
    plt.show()


# Отображение карты вейвлет-коэффициентов
def plot_swt(data, t_d):
    plt.figure()
    coef, freqs = pywt.cwt(data, np.arange(1, 20), 'morl',
                       sampling_period=t_d)
    plt.pcolor(range(1,6401), freqs, coef)
    plt.show()


# Вывод ИХ фильтра по его коэффициентам
def impz(b,a, name):
    impulse = [0]*100
    impulse[0] =1.
    x = range(100)
    response = lfilter(b,a,impulse)
    plt.stem(x, response)
    plt.ylabel('Amplitude')
    plt.xlabel(r'n (samples)')
    plt.title(name)
    plt.show()


# Вывод АЧХ и ФЧХ фильтра по его коэффициентам
def mfreqz(b,a):
    w,h = freqz(b,a)
    h_dB = 20 * np.log10 (abs(h))
    plt.subplot(211)
    plt.plot(w/max(w),h_dB)
    plt.ylabel('Magnitude (db)')
    plt.xlabel(r'Normalized Frequency (x$\pi$rad/sample)')
    plt.title(r'Frequency response')
    plt.grid('on')
    plt.ylim(bottom=-30)
    plt.subplot(212)
    h_Phase = np.unwrap(np.arctan2(np.imag(h),np.real(h)))
    plt.plot(w/max(w),h_Phase)
    plt.ylabel('Phase (radians)')
    plt.xlabel(r'Normalized Frequency (x$\pi$rad/sample)')
    plt.title(r'Phase response')
    plt.grid('on')
    plt.subplots_adjust(hspace=1.5)
    plt.show()


# Вывод АЧХ и ФЧХ фильтра по его коэффициентам
def mfreqz3(b,a, names):
    lines = ["-","--","-.",":"]
    plt.subplot(211)
    for i in range(3):
        w,h = freqz(b[i],a[i])
        h_dB = 20 * np.log10 (abs(h))
        plt.plot(w/max(w),h_dB, linestyle=lines[i])
        plt.legend(names)
    plt.ylabel('Magnitude (db)')
    plt.xlabel(r'Normalized Frequency (x$\pi$rad/sample)')
    plt.title(r'Frequency response')
    plt.grid('on')
    plt.ylim(top=1,bottom=-30)
    plt.subplot(212)
    for i in range(3):
        w,h = freqz(b[i],a[i])
        h_Phase = np.unwrap(np.arctan2(np.imag(h),np.real(h)))
        plt.plot(w/max(w),h_Phase, linestyle=lines[i])
        #plt.legend(names)
    
    plt.ylabel('Phase (radians)')
    plt.xlabel(r'Normalized Frequency (x$\pi$rad/sample)')
    plt.title(r'Phase response')
    plt.grid('on')
    plt.subplots_adjust(hspace=0.5)
    plt.show()


""" Генерация импульсов 
Функции генерации последовательности импульсов и одиночного импульса
"""


# Формирование последовательности чипов
def generate_sequence(sig_type,t_d, n_chips, t_imp, f_low):
#Ф-я формирования набора последовательности сигналов перестраиваемых по
#частоте
# Аргументы:
#   sig_type - тип сигнала
#   t_d - интервал дискретизации
#   n_chips - число отсетов
#   f_low - нижняя частота
#   запустить функцию формирования сигнала в зависимости от заданного типа
    res_sig = np.array([0])
    sig = []
    n_cnts_sig = math.floor(t_imp/t_d)
    if (sig_type == 'chirp'): # тип в/импульс
        f_mod = f_low/5
        for i in range(0, n_chips):
            sig = [0] * n_cnts_sig
            sig[:n_cnts_sig] = get_chirp_pulse(t_d, t_imp, t_imp*0.96, f_low, f_mod)
            res_sig = np.append(res_sig, sig)
    elif (sig_type == 'radio'): # тип р/импульс
        
        # случайная перестройка для р/импульса
        for i in range(0, n_chips):
            sig = [0] * n_cnts_sig
            random_freq = f_low + np.random.randint(5, size=(1))*8/t_imp
            sig[:n_cnts_sig] = get_radio_pulse(t_d, t_imp, t_imp*0.96, random_freq)
            res_sig = np.append(res_sig, sig)
        
    elif (sig_type == 'AM'):
        f_mod = f_low/10
        for i in range(0, n_chips):
            sig = [0] * n_cnts_sig
            random_freq = f_low + np.random.randint(5, size=(1))*(3*f_mod)
            sig[:n_cnts_sig] = get_AM(t_d, t_imp, random_freq, f_mod, 0.5)
            res_sig = np.append(res_sig, sig)
    res_sig = np.delete(res_sig, 0)
    
    return res_sig


# Формирование сигналов (управляющая функция)
def generate_single_chip(signal_type, *signal_data): # запустить функцию формирования сигнала в зависимости от заданного типа
    if (signal_type == 'video'):    # тип в/импульс
        signal = get_video_pulse(*signal_data)
    elif (signal_type == 'radio'):    # тип р/импульс
        signal = get_radio_pulse(*signal_data)
    elif (signal_type == 'AM'):    # тип АМ сигнал  
        signal = get_AM(*signal_data)
    elif (signal_type == 'chirp'):
        signal = get_chirp_pulse(*signal_data)
    return signal


# Функции формирования чипа
# Формирование в/импульса
def get_video_pulse(t_d, t_window, t_imp):
# Ф-я формировнаия в/импульса
# Аргументы:    
#   % t_d - период дискретизации, с t_window - период рассмотрения, с
#   t_imp - длительность импульса, с
    n_signal = math.floor(t_window/t_d)  # число отсчетов в рассматриваемом интерв.
    n_imp = math.floor(t_imp/t_d)        # число отсчетов импульса
    signal = n_signal * [0]            # сформировать пустую выборку
    for i in range(0, n_imp):     # поставить первые n_imp отсчетов
        signal[i] = 1
    return signal


# Формирование р/импульса
def get_radio_pulse(t_d, t_window, t_imp, f_carrier_hz):
# Ф-я формирования р/импульса
# Аргументы:
#   t_d - интервал дискретизации, с    t_window - интервал рассмотрения, с
#   t_imp - длительность импульса, с   f_carrier_hz - частота несущей, Гц
    n_signal = math.floor(t_window/t_d) # число отсчетов в рассматриваемом интерв.
    n_imp = math.floor(t_imp/t_d)       # число отсчетов импульса
    pulse_amp = np.zeros(n_signal)
    for i in range(1, n_imp):     # сформировать огиб.
        pulse_amp[i] = 1
    t_axis = np.linspace(0, t_window, n_signal)   # формирование временной оси
    carrier = np.sin(2*math.pi*f_carrier_hz*t_axis)    # сформировать несущую
    signal = pulse_amp*carrier    # несущая на амплитудные значения
    return signal


# Формирование АМ/сигнала
def get_AM(t_d, t_window, f_carrire_hz, f_mod_hz, i_mod):
# Ф-я формирования АМ сигнала
#   t_d - интервал дискретизации, с    t_window - интервал рассмотрения, с
#   f_carrier_hz - частота несущей, Гц
#   f_mod_hz - частота модуляции, Гц   i_mod - глубина модуляции (Amin/Amax)
    n_signal = math.floor(t_window/t_d)                 # число отсчетов в интервале рассмот.
    t_axis = np.linspace(0, t_window, n_signal)       # ось времени
    am_mult = 1/(2/i_mod) 
    am_shift = 1-am_mult    # множитель для расчета огибающей АМ
    a_modulation = np.sin(2*math.pi*f_mod_hz*t_axis)*am_mult+am_shift  # огибающая
    signal = np.sin(2*math.pi*f_carrire_hz*t_axis) * a_modulation     # формирование АМ сигнала
    return signal


# Формирование ЛЧМ
def get_chirp_pulse(t_d, t_window, t_imp, f_start_hz, f_chirp_hz):
# Ф-я формирования р/импульса
# Аргументы:
#   t_d - интервал дискретизации, с    t_window - интервал рассмотрения, с
#   t_imp - длительность импульса, с   f_center_hz - серединная частота, Гц
#   f_chirp_hz - ширина полосы, Гц
    n_signal = math.floor(t_window/t_d)     # число отсчетов в рассматриваемом интерв.
    n_imp = math.floor(t_imp/t_d)           # число отсчетов импульса
    t_axis = np.linspace(0, t_imp, n_imp)   # формирование временной оси
    chirp_band = f_chirp_hz/t_imp
    chirp_sin_arg = f_start_hz + chirp_band/2*t_axis
    signal = [0] * n_signal
    signal[1:n_imp+1] = np.sin(2*math.pi*chirp_sin_arg*t_axis)     # сформировать несущую
#     signal(n_imp:end) =  sin(2*pi*f_moment(n_imp).*t_axis(n_imp:end))
    return signal


""" Фильтры 
"""


# Применить к сигналу фильтр с "идеальной" АЧХ
def apply_ideal_filter(t_d, filter_type, f_cut_hz, signal_in):
# Ф-я фильтрации фильтром с идеальной АЧХ
# Аргументы:
#   t_d - интервал дискретизации,с 
#   filter_type: 'LP' - НЧ, 'HP' - ВЧ, 'BP' - полосовой, 'S' - заградительный
#   f_cut_hz - частота среза, Гц
#   signal_in - входной сигнал (вектор строка)
    f_d = 1/t_d    # расчет частоты дискретизации
    n_signal = len(signal_in)
    f_axis = np.linspace(0, f_d/2, int(np.floor(n_signal / 2)))
    if (isinstance(f_cut_hz, float)):         # Если одно значение частоты среза, то НЧ или ВЧ фильтр
        f_cut_i_num = list(filter(lambda i: i > f_cut_hz, f_axis))[0]
        f_cut_i = np.argmax(f_axis > f_cut_hz)
        filter_f_char = np.ones(int(np.ceil((n_signal+1)/2)))
        if (filter_type == 'HP'):
            filter_f_char[0:f_cut_i] = 0
            #filter_f_char(1:f_cut_i-1) = 0   # если НЧ 
        elif (filter_type == 'LP'): 
            filter_f_char[f_cut_i:] = 0
            #filter_f_char(f_cut_i:end) = 0   # если ВЧ
        else:   
            print('Ошибка вввода: Тип фильтра')

        # формирование АЧХ для двухстороннего спектра
        # с предпоследнего отсчета
        tmp_filter_f_char = np.flip(filter_f_char[1:])
        if (n_signal % 2 == 0):
            tmp_filter_f_char = tmp_filter_f_char[:len(tmp_filter_f_char) - 1]
        filter_f_char = np.append(filter_f_char, tmp_filter_f_char)
    elif (isinstance(f_cut_hz, list)):     # Если два значения частоты среза, то ПФ или ЗФ фильтр
        #TODO: РЕАЛИЗОВАТЬ ВАРИАНТЫ ПФ И ЗФ
        print("none")
    else:
        print('Ошибка ввода: Граничные частоты фильтра')
        
    signal_in_sp = fft(signal_in, n_signal)
    signal_out_sp = signal_in_sp * filter_f_char
    signal_out = ifft(signal_out_sp, n_signal)
    if (np.prod(np.iscomplex(signal_out))): 
        print('Допущена ошибка при формировании фильтра: сигнал после фильтрации комплексный')
    return signal_out


# Применить к сигналу фильтр Баттерворта
def apply_butt_filter(t_d, filter_type, f_cut_hz, filter_order, signal_in):
# Ф-я фильтрации фильтром с АЧХ Баттерворта
# Аргументы:
#   t_d - интервал дискретизации,с 
#   filter_type: 'LP' - НЧ, 'HP' - ВЧ, 'BP' - полосовой, 'S' - заградительный
#   f_cut_hz - частота среза, Гц
#   filter_order - порядок фильтра
#   signal_in - входной сигнал (вектор строка)
    b, a = create_butt_filter(t_d, filter_type, f_cut_hz, filter_order)  # формирование коэффициентов фильтра Баттерворта
    
    signal_out = lfilter(b,a,signal_in)     # фильтрация сигнала
    return signal_out


# Применить к сигналу фильтр Чебышёва
def apply_cheb2_filter(t_d, filter_type, f_cut_hz, filter_order, bnd_att_db, signal_in):
# Ф-я фильтрации фильтром с АЧХ Чебышева
# Аргументы:
#   t_d - интервал дискретизации,с 
#   filter_type: 'LP' - НЧ, 'HP' - ВЧ, 'BP' - полосовой, 'S' - заградительный
#   f_cut_hz - частота среза, Гц
#   filter_order - порядок фильтра
#   bnd_att_db - подавление боковых лепестков АЧХ
#   signal_in - входной сигнал (вектор строка)
    b, a = create_cheb2_filter(t_d, filter_type, f_cut_hz, filter_order, bnd_att_db)  # формирование фильтра
    
    signal_out = lfilter(b,a,signal_in)                 # фильтрация сигнала
    return signal_out


# Расчет коэффициентов фильтра Баттерворта
def create_butt_filter(t_d, filter_type, f_cut_hz, filter_order):
    # Ф-я расчета коэффициентов передаточной функции фильтра Баттерворта
    # Аргументы:
    #   t_d - интервал дискретизации,с 
    #   filter_type: 'LP' - НЧ, 'HP' - ВЧ, 'BP' - полосовой, 'S' - заградительный
    #   f_cut_hz - частота среза, Гц
    #   filter_order - порядок фильтра
    #   signal_in - входной сигнал (вектор строка)
    if isinstance(f_cut_hz, float):      # если одно значение - ВЧ или НЧ фильтр
        if (filter_type == 'LP'):
            f_name = 'lowpass'
        elif (filter_type == 'HP'):
            f_name = 'highpass'
        else:
            print('Ошибка вввода: Тип фильтра')
    elif (isinstance(f_cut_hz, list)):
        if (filter_type == 'BP'):
            f_name = 'bandpass'
        elif (filter_type == 'S'):
            f_name = 'bandstop'
        else:
            print('Ошибка вввода: Тип фильтра')

    else:
        error('Ошибка ввода: Граничные частоты фильтра')
    
    f_d = 1/t_d            # частота дискретизации
    w_n = f_cut_hz/(f_d/2) # нормированная частота среза
    b, a = butter(filter_order, w_n, f_name)  # формирование коэффициентов фильтра Баттерворта
    return b, a


# Расчет коэффициентов Фильтра Чебышева
def create_cheb2_filter(t_d, filter_type, f_cut_hz, filter_order, bnd_att_db):
# Ф-я расчета коэффициентов передаточной характеристики фильтр с АЧХ Чебышева
# Аргументы:
#   t_d - интервал дискретизации,с 
#   filter_type: 'LP' - НЧ, 'HP' - ВЧ, 'BP' - полосовой, 'S' - заградительный
#   f_cut_hz - частота среза, Гц
#   filter_order - порядок фильтра
#   bnd_att_db - подавление боковых лепестков АЧХ
#   signal_in - входной сигнал (вектор строка)
    if isinstance(f_cut_hz, float):      # если передано одно зачение частоты среза
        if (filter_type == 'LP'):
            f_name = 'lowpass'           # НЧ фильтр
        elif (filter_type == 'HP'):
            f_name = 'highpass'          # ВЧ фильтр
        else:
            print('Ошибка вввода: Тип фильтра') # если передано два значения частот среза
    elif (isinstance(f_cut_hz, list)):
        if (filter_type == 'BP'):
            f_name = 'bandpass'                 # полосовой фильтр
        elif (filter_type == 'S'):
            f_name = 'bandstop'                 # заградительный фильтр
        else:
            print('Ошибка вввода: Тип фильтра')

    else:
        error('Ошибка ввода: Граничные частоты фильтра')

    f_d = 1/t_d            # частота дискретизации, Гц
    w_n = f_cut_hz/(f_d/2) # нормированная частота среза (w/w_d)
    b, a = cheby2(filter_order, bnd_att_db, w_n, f_name)  # формирование фильтра
    return b, a


def create_cheb1_filter(t_d, filter_type, f_cut_hz, filter_order, bnd_att_db):
# Ф-я расчета коэффициентов передаточной характеристики фильтр с АЧХ Чебышева
# Аргументы:
#   t_d - интервал дискретизации,с 
#   filter_type: 'LP' - НЧ, 'HP' - ВЧ, 'BP' - полосовой, 'S' - заградительный
#   f_cut_hz - частота среза, Гц
#   filter_order - порядок фильтра
#   bnd_att_db - подавление боковых лепестков АЧХ
#   signal_in - входной сигнал (вектор строка)
    if isinstance(f_cut_hz, float):      # если передано одно зачение частоты среза
        if (filter_type == 'LP'):
            f_name = 'lowpass'           # НЧ фильтр
        elif (filter_type == 'HP'):
            f_name = 'highpass'          # ВЧ фильтр
        else:
            print('Ошибка вввода: Тип фильтра') # если передано два значения частот среза
    elif (isinstance(f_cut_hz, list)):
        if (filter_type == 'BP'):
            f_name = 'bandpass'                 # полосовой фильтр
        elif (filter_type == 'S'):
            f_name = 'bandstop'                 # заградительный фильтр
        else:
            print('Ошибка вввода: Тип фильтра')

    else:
        error('Ошибка ввода: Граничные частоты фильтра')
        
    f_d = 1/t_d            # частота дискретизации, Гц
    w_n = f_cut_hz/(f_d/2) # нормированная частота среза (w/w_d)
    b, a = cheby1(filter_order, bnd_att_db, w_n, f_name)  # формирование фильтра
    return b, a

