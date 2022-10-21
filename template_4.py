# %%
from toolbox import *

# Задача - выполнить анализ сигнала с использованием спектрограмм и вейвлет-разложения
# Исходные параметры
f_d = 8e6               # частота дискретизации, Гц
t_d = 1/f_d             # период дискретизации
t_impulse = 40e-6       # длительность импульса (t chip)
t_window = 80e-6        # длительность рассматриваемого интервала
f_carrier = 0.80e6      # частота несущей, Гц
f_mod = 0.6e6           
i_mod = 0.8             # параметры модуляции (частота, глубина)
n_chips = 8             
# Формирование последовательности отсчетов
s_c = generate_sequence('chirp', t_d, n_chips, t_impulse, f_carrier);   # ЛЧМ импульс
s_r = generate_sequence('radio', t_d, n_chips, t_impulse, f_carrier);   # р/импульс
s_am = generate_sequence('AM', t_d, n_chips, t_impulse, f_carrier);     # АМ сигнал

# %%
# Анализ последовательности
signal_out = s_r  # рез-татов фильтрации

plot_signal([[t_d, signal_out]])
window_opt_len = len(signal_out)//n_chips
plt.figure()

# Малый размер окна
plt.subplot(311)
powerSpectrum, freqenciesFound, time, imageAxis = plt.specgram(signal_out, 
    NFFT=8,window=np.hamming(8), noverlap=0)
plt.xlabel('Sample')
plt.ylabel('Normalized Frequency')
# Большой размер окна с шагом на один отсчёт
plt.subplot(312)
powerSpectrum, freqenciesFound, time, imageAxis = plt.specgram(signal_out, 
    NFFT=window_opt_len, window=np.hamming(window_opt_len), noverlap=window_opt_len-1)
plt.xlabel('Sample')
plt.ylabel('Normalized Frequency')
# Большой размер окна с шагом на величину окна
plt.subplot(313)
powerSpectrum, freqenciesFound, time, imageAxis = plt.specgram(signal_out, 
    NFFT=window_opt_len, window=np.hamming(window_opt_len), noverlap=0)
plt.xlabel('Sample')
plt.ylabel('Normalized Frequency')
plt.show()

plt.figure()
cwtmatr, freqs = pywt.cwt(signal_out, np.arange(1,63), 'morl')
plt.imshow(cwtmatr, extent=[-1, 1, 1, 63], cmap='PRGn', aspect='auto',
           vmax=abs(cwtmatr).max(), vmin=-abs(cwtmatr).max())  
plt.show() 

