# %%
#!/usr/bin/env python
# coding: utf-8
import matplotlib.pyplot as plt
import math
import numpy as np
from scipy.fftpack import fft, fftfreq, fftshift, ifft, rfft
from scipy.signal import butter, lfilter, cheby2, cheby1, impulse2, freqz


# Служебные функции общего назначения

# Формирование сигналов (управляющая функция)
def generate_single_chirp(signal_type, *signal_data): # запустить функцию формирования сигнала в зависимости от заданного типа
    if (signal_type == 'video'):    # тип в/импульс
        signal = get_video_pulse(*signal_data)
    elif (signal_type == 'radio'):    # тип р/импульс
        signal = get_radio_pulse(*signal_data)
    elif (signal_type == 'AM'):    # тип АМ сигнал  
        signal = get_AM(*signal_data)
    return signal


# Отображение сигналов
def plot_signal(signal_args):
#PLOT_SIGNAL Отобразить сигнал, в случае нескольких сигналов, они
#передаются по горизонтали
# Аргументы:
#   один аргумент - значения сигнала, тогда по оси абсцис отсчеты
#   два аргумента - {временная шкала или шаг дискретизации, значения амплитуд}
#   три аргумента - третьим передается название окна
    lines = ["-","--","-.",":"]
    leg = [];
     
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
            signal = signal_args[i][1];
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
    plt.figure(figsize=(15, 10))
    plt.show()


# Отображение спектр сигнала
def plot_spectum(signal_array):
    t_d_us = signal_array[0][0]*1e6;
    
    fig, axs = plt.subplots(len(signal_array), 1)
    
    for i in range(len(signal_array)):
        signal = signal_array[i][1]
    #[N_signals, n_signal] = size(signal);
        n_signal = signal.size;
        t_window = n_signal*t_d_us;
        f_step_mhz = 1/t_window;
        if n_signal % 2 == 0:
            f_axis = np.arange(0,(n_signal+1)*f_step_mhz/2,f_step_mhz)
        else:
            print('нечетное число отсчетов');
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
            

# Функции формирования чипа

# Формирование в/импульса
def get_video_pulse(t_d, t_window, t_imp):
# Ф-я формировнаия в/импульса
# Аргументы:    
#   % t_d - период дискретизации, с; t_window - период рассмотрения, с
#   t_imp - длительность импульса, с
    n_signal = math.floor(t_window/t_d)  # число отсчетов в рассматриваемом интерв.
    n_imp = math.floor(t_imp/t_d)        # число отсчетов импульса
    signal = n_signal * [0]            # сформировать пустую выборку
    for i in range(0, n_imp):     # поставить первые n_imp отсчетов
        signal[i] = 1;
    return signal


# Формирование р/импульса
def get_radio_pulse(t_d, t_window, t_imp, f_carrier_hz):
# Ф-я формирования р/импульса
# Аргументы:
#   t_d - интервал дискретизации, с;    t_window - интервал рассмотрения, с
#   t_imp - длительность импульса, с;   f_carrier_hz - частота несущей, Гц
    n_signal = math.floor(t_window/t_d); # число отсчетов в рассматриваемом интерв.
    n_imp = math.floor(t_imp/t_d);       # число отсчетов импульса
    pulse_amp = np.zeros(n_signal);
    for i in range(1, n_imp):     # сформировать огиб.
        pulse_amp[i] = 1;
    t_axis = np.linspace(0, t_window, n_signal);   # формирование временной оси
    carrier = np.sin(2*math.pi*f_carrier_hz*t_axis);    # сформировать несущую
    signal = pulse_amp*carrier;    # несущая на амплитудные значения
    return signal


# Формирование АМ/сигнала
def get_AM(t_d, t_window, f_carrire_hz, f_mod_hz, i_mod):
# Ф-я формирования АМ сигнала
#   t_d - интервал дискретизации, с;    t_window - интервал рассмотрения, с
#   f_carrier_hz - частота несущей, Гц;
#   f_mod_hz - частота модуляции, Гц;   i_mod - глубина модуляции (Amin/Amax)
    n_signal = math.floor(t_window/t_d);                 # число отсчетов в интервале рассмот.
    t_axis = np.linspace(0, t_window, n_signal);       # ось времени
    am_mult = 1/(2/i_mod); 
    am_shift = 1-am_mult;    # множитель для расчета огибающей АМ
    a_modulation = np.sin(2*math.pi*f_mod_hz*t_axis)*am_mult+am_shift;  # огибающая
    signal = np.sin(2*math.pi*f_carrire_hz*t_axis) * a_modulation;     # формирование АМ сигнала
    return signal


# Функции формирования фильтров

# Применить к сигналу фильтр с "идеальной" АЧХ
def apply_ideal_filter(t_d, filter_type, f_cut_hz, signal_in):
# Ф-я фильтрации фильтром с идеальной АЧХ
# Аргументы:
#   t_d - интервал дискретизации,с; 
#   filter_type: 'LP' - НЧ, 'HP' - ВЧ, 'BP' - полосовой, 'S' - заградительный
#   f_cut_hz - частота среза, Гц;
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
        elif (filter_type == 'LP'): 
            filter_f_char[f_cut_i:] = 0
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
#   t_d - интервал дискретизации,с; 
#   filter_type: 'LP' - НЧ, 'HP' - ВЧ, 'BP' - полосовой, 'S' - заградительный
#   f_cut_hz - частота среза, Гц;
#   filter_order - порядок фильтра
#   signal_in - входной сигнал (вектор строка)
    b, a = create_butt_filter(t_d, filter_type, f_cut_hz, filter_order);  # формирование коэффициентов фильтра Баттерворта
    
    signal_out = lfilter(b,a,signal_in);     # фильтрация сигнала
    return signal_out


# Применить к сигналу фильтр Чебышёва
def apply_cheb2_filter(t_d, filter_type, f_cut_hz, filter_order, bnd_att_db, signal_in):
# Ф-я фильтрации фильтром с АЧХ Чебышева
# Аргументы:
#   t_d - интервал дискретизации,с; 
#   filter_type: 'LP' - НЧ, 'HP' - ВЧ, 'BP' - полосовой, 'S' - заградительный
#   f_cut_hz - частота среза, Гц;
#   filter_order - порядок фильтра;
#   bnd_att_db - подавление боковых лепестков АЧХ;
#   signal_in - входной сигнал (вектор строка)
    b, a = create_cheb2_filter(t_d, filter_type, f_cut_hz, filter_order, bnd_att_db);  # формирование фильтра
    
    signal_out = lfilter(b,a,signal_in);                 # фильтрация сигнала
    return signal_out


# Расчет коэффициентов фильтра Баттерворта
def create_butt_filter(t_d, filter_type, f_cut_hz, filter_order):
    # Ф-я расчета коэффициентов передаточной функции фильтра Баттерворта
    # Аргументы:
    #   t_d - интервал дискретизации,с; 
    #   filter_type: 'LP' - НЧ, 'HP' - ВЧ, 'BP' - полосовой, 'S' - заградительный
    #   f_cut_hz - частота среза, Гц;
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
        print('Ошибка ввода: Граничные частоты фильтра');
    
    f_d = 1/t_d;            # частота дискретизации
    w_n = f_cut_hz/(f_d/2); # нормированная частота среза
    b, a = butter(filter_order, w_n, f_name);  # формирование коэффициентов фильтра Баттерворта
    return b, a


# Расчет коэффициентов Фильтра Чебышева
def create_cheb2_filter(t_d, filter_type, f_cut_hz, filter_order, bnd_att_db):
# Ф-я расчета коэффициентов передаточной характеристики фильтр с АЧХ Чебышева
# Аргументы:
#   t_d - интервал дискретизации,с; 
#   filter_type: 'LP' - НЧ, 'HP' - ВЧ, 'BP' - полосовой, 'S' - заградительный
#   f_cut_hz - частота среза, Гц;
#   filter_order - порядок фильтра;
#   bnd_att_db - подавление боковых лепестков АЧХ;
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
        print('Ошибка ввода: Граничные частоты фильтра')

    f_d = 1/t_d;            # частота дискретизации, Гц
    w_n = f_cut_hz/(f_d/2); # нормированная частота среза (w/w_d)
    b, a = cheby2(filter_order, bnd_att_db, w_n, f_name);  # формирование фильтра
    return b, a


def create_cheb1_filter(t_d, filter_type, f_cut_hz, filter_order, bnd_att_db):
# Ф-я расчета коэффициентов передаточной характеристики фильтр с АЧХ Чебышева
# Аргументы:
#   t_d - интервал дискретизации,с; 
#   filter_type: 'LP' - НЧ, 'HP' - ВЧ, 'BP' - полосовой, 'S' - заградительный
#   f_cut_hz - частота среза, Гц;
#   filter_order - порядок фильтра;
#   bnd_att_db - подавление боковых лепестков АЧХ;
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
        print('Ошибка ввода: Граничные частоты фильтра')
        
    f_d = 1/t_d;            # частота дискретизации, Гц
    w_n = f_cut_hz/(f_d/2); # нормированная частота среза (w/w_d)
    b, a = cheby1(filter_order, bnd_att_db, w_n, f_name);  # формирование фильтра
    return b, a


# %%
# Лабораторная работа №3
# Задача - сформировать сигналы и выполнить их фильтрацию
# Дополнить код реализациями полосового и заградительного фильтров
# Формирование сигнала
f_d = 1e9              # частота дискретизации, Гц
t_d = 1/f_d            # период дискретизации
t_win = 2e-7           # длительность рассматриваемого интервала
t_ch = 8e-9            # длительность импульса (t chirp)
f_car = 140e6          # частота несущей, Гц
f_mod = 5e6 
i_mod = 0.5   # параметры модуляции


# Формирование сигналов
s_v = generate_single_chirp('video', t_d, t_win, t_ch)                 # в/импульс
s_r = generate_single_chirp('radio', t_d, t_win, t_ch, f_car)          # р/импульс
s_am = generate_single_chirp('AM', t_d, t_win, f_car, f_mod, i_mod)    # АМ сигнал


# Отображение данных
signal = np.roll(s_v, math.floor(len(s_r)/2)) # выбрать сигнал для отображения
  
plot_signal([[t_d,signal]])     # отобразить сигнал
plot_spectum([[t_d, signal]]);   # отобразить спектр сигнала


# Сравнение фильтров
# идеальный КИХ фильтр
signal_in = signal     # выбрать сигнал для дальнейшей обработки
f_type = 'LP'          # тип фильтра (LP/HP/BP/S)
f_cut_hz = 1/t_ch      # частота среза, Гц
but_order = 3  
cheb_order = 5 # порядок фильтров
s_band_attenuation_db = 20     # внеполосное ослабление лепестков для фильтра Чебышева


signal_out_idl = apply_ideal_filter(t_d, f_type, f_cut_hz, signal_in);              # идеальный ф.
signal_out_btr = apply_butt_filter(t_d, f_type, f_cut_hz, but_order, signal_in);    # ф. Баттерворта
signal_out_chb = apply_cheb2_filter(t_d, f_type, f_cut_hz, cheb_order, s_band_attenuation_db, signal_in);              # ф. Чебышева


# Отображение данных
signal_out = [[t_d, signal_out_idl, 'ideal'], [t_d, signal_out_btr, 'butterworth'], 
              [t_d, signal_out_chb, 'chebyshev']]  # конкатенация рез-татов фильтрации
plot_signal(signal_out)                 # отобразить фильтрованные сигналы
plot_spectum(signal_out)       # отобразить спектр (ф-я возвращает указатель на объект)


def impz(b,a):
    impulse = [0]*100
    impulse[0] =1.
    x = range(100)
    response = lfilter(b,a,impulse)
    plt.stem(x, response)
    plt.ylabel('Amplitude')
    plt.xlabel(r'n (samples)')
    plt.title(r'Impulse response')


def mfreqz(b,a):
    w,h = freqz(b,a)
    h_dB = 20 * np.log10 (abs(h))
    plt.subplot(211)
    plt.plot(w/max(w),h_dB)
    plt.ylabel('Magnitude (db)')
    plt.xlabel(r'Normalized Frequency (x$\pi$rad/sample)')
    plt.title(r'Frequency response')
    plt.grid('on')
    plt.subplot(212)
    h_Phase = np.unwrap(np.arctan2(np.imag(h),np.real(h)))
    plt.plot(w/max(w),h_Phase)
    plt.ylabel('Phase (radians)')
    plt.xlabel(r'Normalized Frequency (x$\pi$rad/sample)')
    plt.title(r'Phase response')
    plt.grid('on')
    plt.subplots_adjust(hspace=0.5)
    plt.show()


# Вывод АЧХ и ФЧХ
b, a = create_cheb2_filter(t_d, f_type, f_cut_hz, cheb_order, s_band_attenuation_db);      # формирование фильтра
n_imp_resp = 2^7;   # длительность интервала (в отсчетах) для которого строится имп. характеристика
n_freq_resp = 2^10;  # число отсчетов АЧХ и ФЧХ
   
impz(b,a)

mfreqz(b, a)







# %%
