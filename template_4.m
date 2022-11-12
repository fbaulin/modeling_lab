% Исходные параметры
f_d = 20e6;                          % частота дискретизации, Гц
t_d = 1/f_d;                        % период дискретизации
t_impulse = 20e-6;                % длительность импульса (t chip)
t_window = 80e-6;                % длительность рассматриваемого интервала
f_carrier = 1.80e6;                  % частота несущей, Гц
f_mod = 0.6e6; i_mod = 0.8; % параметры модуляции (частота, глубина)
%% Работа №3 Моделирование фильтров
% Задача - выполнить анализ сигнала с использованием спектрограмм и вейвлет-разложения
% Формирование одиночных импульсов
% s_c = generate_single_chip('chirp', t_d, t_window, t_impulse, f_carrier, f_mod);    % ЛЧМ импульс
s_v = generate_single_chip('video', t_d, t_window, t_impulse);                      % в/импульс
% s_r = generate_single_chip('radio', t_d, t_window, t_impulse, f_carrier);           % р/импульс
% s_am = generate_single_chip('AM', t_d, t_window, f_carrier, f_mod, i_mod);          % АМ сигнал
%% Отображение импульса
signal = circshift(s_v,floor(length(s_v)/2));  % выбрать сигнал для отображения
plot_signal(t_d,signal)     % отобразить сигнал
plot_spectum(t_d,signal);   % отобразить спектр сигнала
%% Сравнение фильтров
% Исходные параметры
signal_in = signal;         % выбрать сигнал для дальнейшей обработки
f_type = 'HP';              % тип фильтра (LP/HP/BP/S)
f_cut_hz = 1/t_impulse;             % частота среза, Гц
but_order = 3;  cheb_order = 3;     % порядок фильтров
s_band_attenuation_db = 20;         % внеполосное ослабление лепестков для фильтра Чебышева

signal_out_idl = apply_ideal_filter(t_d, f_type, f_cut_hz, signal_in);              % идеальный ф.
signal_out_btr = apply_butt_filter(t_d, f_type, f_cut_hz, but_order, signal_in);    % ф. Баттерворта
signal_out_chb = apply_cheb2_filter(t_d, f_type, f_cut_hz, cheb_order, ...
    s_band_attenuation_db, signal_in);                                              % ф. Чебышева
% Отображение результатов
signal_out = [signal_out_idl; signal_out_btr; signal_out_chb];  % конкатенация рез-татов фильтрации
plot_signal(t_d,signal_out)                 % отобразить фильтрованные сигналы
legend('ideal','butterworth','chebyshev')   % подписать графики
hndls = plot_spectum(t_d,signal_out);       % отобразить спектр (ф-я возвращает указатель на объект)
% подписать графики
legend(hndls(1),'ideal');
legend(hndls(2),'butterworth');
legend(hndls(3),'chebyshev');
%% Работа №4 Частотно-временной анализ
% Формирование последовательности отсчетов
n_chips = 8;
s_c = generate_sequence('chirp', t_d, n_chips, t_impulse, f_carrier, f_mod);   % ЛЧМ импульс
s_r = generate_sequence('radio', t_d, n_chips, t_impulse, f_carrier);   % р/импульс
s_am = generate_sequence('AM', t_d, n_chips, t_impulse, f_carrier, f_mod);     % АМ сигнал
plot_signal(t_d,[s_c;s_r;s_am])

%% Отображение последовательности ЛЧМ импульсов
figure('Name', 'СПГ ЛЧМ','WindowStyle','docked','Color','white')
subplot(3,1,1); spectrogram(s_c,6,2,'yaxis');       % спектрограмма (короткое окно)
subplot(3,1,2); spectrogram(s_c,100,80,'yaxis');    % спектрограмма (длинное окно, короткий шаг)
subplot(3,1,3); spectrogram(s_c,100,0,'yaxis');     % спектрограмма (длинное окно)
figure('Name', 'Вейв. ЛЧМ','WindowStyle','docked','Color','white')
cwt(s_c,f_d,'amor');

%% Отображение последовательности радиоимпульсов
figure('Name', 'СПГ Р.имп.','WindowStyle','docked','Color','white')
subplot(3,1,1); spectrogram(s_r,6,2,n_fft,'yaxis');     % спектрограмма (короткое окно)
subplot(3,1,2); spectrogram(s_r,400,390,400,'yaxis');   % спектрограмма (длинное окно, короткий шаг)
subplot(3,1,3); spectrogram(s_r,400,0,400,'yaxis');     % спектрограмма (длинное окно)
figure('Name', 'Вейв. Р.имп','WindowStyle','docked','Color','white')
cwt(s_r,f_d,'amor');

%% Отображение последовательности АМ чипов
figure('Name', 'СПГ АМ чип','WindowStyle','docked','Color','white')
subplot(3,1,1); spectrogram(s_am,6,2,'yaxis');      % спектрограмма (короткое окно)
subplot(3,1,2); spectrogram(s_am,400,80,'yaxis');   % спектрограмма (длинное окно, короткий шаг)
subplot(3,1,3); spectrogram(s_am,400,0,'yaxis');    % спектрограмма (длинное окно)
figure('Name', 'Вейв. АМ чип','WindowStyle','docked','Color','white')
cwt(s_am,f_d);

%% Служебные функции общего назначение
% Формирование сигналов (управляющая функция)
function signal = generate_single_chip(type,varargin)
    switch type     % запустить функцию формирования сигнала в зависимости от заданного типа
        case 'video'    % тип в/импульс
            signal = get_video_pulse(varargin{:});
        case 'radio'    % тип р/импульс
            signal = get_radio_pulse(varargin{:});
        case 'chirp'
            signal = get_chirp_pulse(varargin{:});
        case 'AM'       % тип АМ сигнал
            signal = get_AM(varargin{:});
    end
end
% Формирование последовательности чипов
function signal = generate_sequence(type,t_d, n_chips, t_imp, f_low, varargin)
%Ф-я формирования набора последовательности сигналов перестраиваемых по
%частоте
% Аргументы:
%   type - тип сигнала
%   t_d - интервал дискретизации
%   n_chips - число импульсов
%   t_imp - длительность одного импульса
%   f_low - нижняя частота
% Дополнительные аргументы:
%   zастота модуляции
    switch type     % запустить функцию формирования сигнала в зависимости от заданного типа
        case 'chirp'    % тип в/импульс
            f_mod =  varargin{1};
            get_chip = @(i) get_chirp_pulse(t_d, t_imp, t_imp, f_low, f_mod);
        case 'radio'    % тип р/импульс
            random_freq = f_low + randi(5,1,n_chips).*8/t_imp;       % случайная перестройка для р/импульса
            get_chip = @(i) get_radio_pulse(t_d, t_imp, t_imp, random_freq(i));
        case 'AM'       % тип АМ сигнал
            f_mod = varargin{1};   % приравнять частоту модуляции 1/10 частоте несущей
            random_freq = f_low + (3*f_mod)*randi(5,1,n_chips);    % формирование случайных частот
            get_chip = @(i) get_AM(t_d, t_imp, random_freq(i), f_mod, 0.5); 
    end
    n_cnts_signal = floor(t_imp/t_d);
    signal = zeros(n_cnts_signal, n_chips);
    
    for i_chip = 1:n_chips
        signal(:,i_chip)=get_chip(i_chip);
    end
    signal(end-8:end,:)=0;
    signal = signal(:).';
end
% Отображение сигналов
function [varargout] = plot_signal(varargin)
%PLOT_SIGNAL Отобразить сигнал, в случае нескольких сигналов, они
%передаются по горизонтали
% Аргументы:
%   один аргумент - значения сигнала, тогда по оси абсцис отсчеты
%   два аргумента - {временная шкала или шаг дискретизации, значения амплитуд}
    switch nargin
        case 1          % если передан только сигнал
            signal = varargin{1};
            x_axis = 0:length(signal)-1;    % временная ось в отсчетах
        case 2          % если переданы временные параметры, сигнал и 
            signal = varargin{2};
            if length(varargin{1})==1       % если передан интервал дискретизации
                x_axis = 0:varargin{1}:varargin{1}*(length(signal)-1);
            else;   x_axis = varargin{1};   % если передана ось времени
            end
        otherwise
            error('Неверный формат ввода: число аргументов');
    end
    h = zeros(size(signal,1));              % вектор с указателями на графики
    markers = {'-r','--g',':b','-.k','-c','--m',':y','-.b'};    % настройки линий
    figure('Name','Сигнал', 'WindowStyle','docked', 'Color','white'); hold on;
    for i=1:size(signal,1)
        h(i) = plot(x_axis*1e6,signal(i,:),markers{i});         % график с сохран. указателя
    end
    xlabel('t, мкc');
    if nargout==1; varargout{1}=h; end  % если предполагается вывод
end
% Отображение спектр сигнала
function [varargout] = plot_spectum(t_d,signal)
% Отображение спектров в одном окне на графиках друг под другом
% Аргументы:
%   t_d - интервал дискретизации
%   signal - отсчеты сигнала (действительного)
    t_d_us = t_d*1e6;
    [N_signals, n_signal] = size(signal);
    n_fft = 10*n_signal;
    t_window_us = n_signal*t_d_us;
    f_step_mhz = 1/t_window_us*n_signal/n_fft;
    if ~mod(n_signal,2)
        f_axis = 0:f_step_mhz:n_fft*f_step_mhz/2;
        f_axis = [-f_axis(end-1:-1:2), f_axis];
    else
        warning('Нечетное число отсчетов');
        f_axis = 0:f_step_mhz:(n_fft-1)*f_step_mhz/2;
        f_axis = [-f_axis(end:-1:2), f_axis];
    end
    signal_spectrum = 1/n_signal*fft(signal,n_fft,2);
    h = zeros(1,N_signals);
    figure('Name','Cпектр', 'WindowStyle','docked', 'Color','white')
    title(['Частота дискретизации ' num2str(1/t_d_us) ' МГц'])
    for i=1:N_signals
        subplot(N_signals, 1, i)
        h(i) = plot(f_axis,fftshift(abs(signal_spectrum(i,:))));
        xlabel('Частота, МГц')
    end  
    if nargout==1; varargout{1} = h; end
end

%% Функции формирования чипа
% Формирование в/импульса
function signal = get_video_pulse(t_d, t_window, t_imp)
% Ф-я формировнаия в/импульса
% Аргументы:    
%   % t_d - период дискретизации, с; t_window - период рассмотрения, с
%   t_imp - длительность импульса, с
    n_signal = floor(t_window/t_d); % число отсчетов в рассматриваемом интерв.
    n_imp = floor(t_imp/t_d);       % число отсчетов импульса
    signal = zeros(1,n_signal);     % сформировать пустую выборку
    signal(1:n_imp)=1;              % поставить первые n_imp отсчетов
end
% Формирование р/импульса
function signal = get_radio_pulse(t_d, t_window, t_imp, f_carrier_hz)
% Ф-я формирования р/импульса
% Аргументы:
%   t_d - интервал дискретизации, с;    t_window - интервал рассмотрения, с
%   t_imp - длительность импульса, с;   f_carrier_hz - частота несущей, Гц
    n_signal = floor(t_window/t_d); % число отсчетов в рассматриваемом интерв.
    n_imp = floor(t_imp/t_d);       % число отсчетов импульса
    pulse_amp = zeros(1, n_signal);     pulse_amp(1:n_imp)=1; % сформировать огиб.
    t_axis = linspace(0, t_window, n_signal);   % формирование временной оси
    carrier = sin(2*pi*f_carrier_hz*t_axis);    % сформировать несущую
    signal = pulse_amp.*carrier;    % несущая на амплитудные значения
end
% Формирование АМ/сигнала
function signal = get_AM(t_d, t_window, f_carrire_hz, f_mod_hz, i_mod)
% Ф-я формирования АМ сигнала
%   t_d - интервал дискретизации, с;    t_window - интервал рассмотрения, с
%   f_carrier_hz - частота несущей, Гц;
%   f_mod_hz - частота модуляции, Гц;   i_mod - глубина модуляции (Amin/Amax)
    n_signal = floor(t_window/t_d);                 % число отсчетов в интервале рассмот.
    t_axis = linspace(0, t_window, n_signal);       % ось времени
    am_mult = 1/(2/i_mod); am_shift = 1-am_mult;    % множитель для расчета огибающей АМ
    a_modulation = sin(2*pi*f_mod_hz*t_axis)*am_mult+am_shift;  % огибающая
    signal = sin(2*pi*f_carrire_hz*t_axis) .* a_modulation;     % формирование АМ сигнала
end
% Формирование ЛЧМ
function signal = get_chirp_pulse(t_d, t_window, t_imp, f_start_hz, f_chirp_hz)
% Ф-я формирования р/импульса
% Аргументы:
%   t_d - интервал дискретизации, с;    t_window - интервал рассмотрения, с
%   t_imp - длительность импульса, с;   f_center_hz - серединная частота, Гц
%   f_chirp_hz - ширина полосы, Гц
    n_signal = floor(t_window/t_d);     % число отсчетов в рассматриваемом интерв.
    n_imp = floor(t_imp/t_d);           % число отсчетов импульса
    t_axis = linspace(0, t_imp, n_imp);   % формирование временной оси
    chirp_band = f_chirp_hz/t_imp;
    chirp_sin_arg = f_start_hz + chirp_band/2*t_axis;
    signal = zeros(1,n_signal);
    signal(1:n_imp) = sin(2*pi*chirp_sin_arg.*t_axis);     % сформировать несущую
%     signal(n_imp:end) =  sin(2*pi*f_moment(n_imp).*t_axis(n_imp:end));     
    
end

%% Функции формирования фильтров
% Формирование фильтра с идеальной АЧХ и заданной(-ыми) частотами среза
function signal_out = apply_ideal_filter(t_d, filter_type, f_cut_hz, signal_in)
% Ф-я фильтрации фильтром с идеальной АЧХ
% Аргументы:
%   t_d - интервал дискретизации,с; 
%   filter_type: 'LP' - НЧ, 'HP' - ВЧ, 'BP' - полосовой, 'S' - заградительный
%   f_cut_hz - частота среза, Гц;
%   signal_in - входной сигнал (вектор строка)
    f_d = 1/t_d;    % расчет частоты дискретизации
    n_signal = length(signal_in);
    f_axis = linspace(0,f_d/2,n_signal/2);
    if length(f_cut_hz)==1         % Если одно значение частоты среза, то НЧ или ВЧ фильтр
        [~, f_cut_i] = find(f_axis>=f_cut_hz,1);        % определение правой границы
        filter_f_char = ones(1,ceil((n_signal+1)/2));   % для 2n -> n+1; для 2n+1 -> n+1
        if      strcmp(filter_type,'HP'); filter_f_char(1:f_cut_i-1) = 0;   % если НЧ 
        elseif  strcmp(filter_type,'LP'); filter_f_char(f_cut_i:end) = 0;   % если ВЧ
        else;   error('Ошибка вввода: Тип фильтра'); 
        end
        % формирование АЧХ для двухстороннего спектра
        filter_f_char = [...
            filter_f_char, ...
            filter_f_char(end-(mod(n_signal,2)==0):-1:2)]; % с предпоследнего отсчета
    elseif length(f_cut_hz)==2     % Если два значения частоты среза, то ПФ или ЗФ фильтр
        %TODO: РЕАЛИЗОВАТЬ ВАРИАНТЫ ПФ И ЗФ
    else; error('Ошибка ввода: Граничные частоты фильтра')
    end
    signal_in_sp = fft(signal_in, n_signal, 2);     % БПФ сигнала
    signal_out_sp = signal_in_sp .* filter_f_char;  % фильтрация
    signal_out = ifft(signal_out_sp, n_signal, 2);  % обратное БПФ спектра сигнала
    % проверка выходного сигнала: если он комплексный, то допущена ошибка
    if ~isreal(signal_out)  
        warning('Допущена ошибка при формировании фильтра: сигнал после фильтрации комплексный')
    end
end
% Формирование фильтра Баттерворта с заднной частотой среза и порядком
function signal_out = apply_butt_filter(t_d, filter_type, f_cut_hz, filter_order, signal_in)
% Ф-я фильтрации фильтром с АЧХ Баттерворта
% Аргументы:
%   t_d - интервал дискретизации,с; 
%   filter_type: 'LP' - НЧ, 'HP' - ВЧ, 'BP' - полосовой, 'S' - заградительный
%   f_cut_hz - частота среза, Гц;
%   filter_order - порядок фильтра
%   signal_in - входной сигнал (вектор строка)
    if length(f_cut_hz)==1      % если одно значение - ВЧ или НЧ фильтр
        switch filter_type
            case 'LP'; f_name = 'low';  % фильтр НЧ
            case 'HP'; f_name = 'high'; % фильтр ВЧ
            otherwise; error('Ошибка вввода: Тип фильтра'); 
        end
    elseif length(f_cut_hz)==2  % если два значения частоты - ПФ ил ЗФ
        switch filter_type
            case 'BP'; f_name = 'bandpass';     % полосовой фильтр
            case 'S'; f_name = 'stop';          % заградительный фильтр
            otherwise; error('Ошибка вввода: Тип фильтра'); 
        end
    else; error('Ошибка ввода: Граничные частоты фильтра');
    end
    f_d = 1/t_d;            % частота дискретизации
    w_n = f_cut_hz/(f_d/2); % нормированная частота среза
    [b,a] = butter(filter_order, w_n, f_name);  % формирование коэффициентов фильтра Баттерворта
    signal_out = filter(b,a,signal_in.').';     % фильтрация сигнала
end
% Формирование фильтра чебышева с заданной частотой среза и порядком
% (внеполосные колебания)
function signal_out = apply_cheb2_filter(t_d, filter_type, f_cut_hz, filter_order, bnd_att_db, signal_in)
% Ф-я фильтрации фильтром с АЧХ Баттерворта
% Аргументы:
%   t_d - интервал дискретизации,с; 
%   filter_type: 'LP' - НЧ, 'HP' - ВЧ, 'BP' - полосовой, 'S' - заградительный
%   f_cut_hz - частота среза, Гц;
%   filter_order - порядок фильтра;
%   bnd_att_db - подавление боковых лепестков АЧХ;
%   signal_in - входной сигнал (вектор строка)
    if length(f_cut_hz)==1  % если передано одно зачение частоты среза
        switch filter_type
            case 'LP'; f_name = 'low';  % НЧ фильтр
            case 'HP'; f_name = 'high'; % ВЧ фильтр
            otherwise; error('Ошибка вввода: Тип фильтра'); 
        end
    elseif length(f_cut_hz)==2  % если передано два значения частот среза
        switch filter_type
            case 'BP'; f_name = 'bandpass'; % полосовой фильтр
            case 'S'; f_name = 'stop';      % заградительный фильтр
            otherwise; error('Ошибка вввода: Тип фильтра'); 
        end
    else; error('Ошибка ввода: Граничные частоты фильтра');
    end
    f_d = 1/t_d;            % частота дискретизации, Гц
    w_n = f_cut_hz/(f_d/2); % нормированная частота среза (w/w_d)
    [b,a] = cheby2(filter_order, bnd_att_db, w_n, f_name);  % формирование фильтра
    signal_out = filter(b,a,signal_in.').';                 % фильтрация сигнала
end
