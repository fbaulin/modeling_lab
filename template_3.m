%% Лабораторная работа №3
% Задача - сформировать сигналы и выполнить их фильтрацию
% Дополнить код реализациями полосового и заградительного фильтров
% Формирование сигнала
f_d = 1e9;              % частота дискретизации, Гц
t_d = 1/f_d;            % период дискретизации
t_win = 2e-7;           % длительность рассматриваемого интервала
t_ch = 8e-9;            % длительность импульса (t chirp)
f_car = 140e6;          % частота несущей, Гц
f_mod = 5e6; i_mod = 0.5;   % параметры модуляции (частота,
% Формирование сигналов
s_v = generate_single_chip('video', t_d, t_win, t_ch);                 % в/импульс
s_r = generate_single_chip('radio', t_d, t_win, t_ch, f_car);          % р/импульс
s_am = generate_single_chip('AM', t_d, t_win, f_car, f_mod, i_mod);    % АМ сигнал
% Отображение данных
signal = circshift(s_v,floor(length(s_r)/2));   % выбрать сигнал для отображения
plot_signal(t_d,signal)     % отобразить сигнал
plot_spectum(t_d,signal);   % отобразить спектр сигнала
%% Сравнение фильтров
% идеальный КИХ фильтр
signal_in = signal;     % выбрать сигнал для дальнейшей обработки
f_type = 'LP';          % тип фильтра (LP/HP/BP/S)
f_cut_hz = 1/t_ch;      % частота среза, Гц
but_order = 3;  cheb_order = 5; % порядок фильтров
s_band_attenuation_db = 20;     % внеполосное ослабление лепестков для фильтра Чебышева

signal_out_idl = apply_ideal_filter(t_d, f_type, f_cut_hz, signal_in);              % идеальный ф.
signal_out_btr = apply_butt_filter(t_d, f_type, f_cut_hz, but_order, signal_in);    % ф. Баттерворта
signal_out_chb = apply_cheb2_filter(t_d, f_type, f_cut_hz, ...
    cheb_order, s_band_attenuation_db, signal_in);              % ф. Чебышева
% Отображение данных
signal_out = [signal_out_idl; signal_out_btr; signal_out_chb];  % конкатенация рез-татов фильтрации
plot_signal(t_d,signal_out)                 % отобразить фильтрованные сигналы
legend('ideal','butterworth','chebyshev')   % подписать графики
hndls = plot_spectum(t_d,signal_out);       % отобразить спектр (ф-я возвращает указатель на объект)
legend(hndls(1),'ideal');  legend(hndls(2),'butterworth');  legend(hndls(3),'chebyshev');   % подписать графики
%% Вывод АЧХ и ФЧХ
[b, a] = create_cheb2_filter(t_d, f_type, f_cut_hz, cheb_order, s_band_attenuation_db);      % формирование фильтра
n_imp_resp = 2^7;   % длительность интервала (в отсчетах) для которого строится имп. характеристика
n_freq_resp = 2^10;  % число отсчетов АЧХ и ФЧХ
figure('Name', 'ИХ','WindowStyle','docked','Color','white');   
impz(b, a, n_imp_resp, f_d)     % построение ИХ (можно вытащить отсчеты - см. справку)
figure('Name', 'АЧХ,ФЧХ','WindowStyle','docked','Color','white');   
freqz(b, a, n_freq_resp, f_d)   % построение АЧХ и ФЧХ (можно вытащить отсчеты - см. справку)

%% Служебные функции общего назначение
% Формирование сигналов (управляющая функция)
function signal = generate_single_chip(type,varargin)
    switch type     % запустить функцию формирования сигнала в зависимости от заданного типа
        case 'video'    % тип в/импульс
            signal = get_video_pulse(varargin{:});
        case 'radio'    % тип р/импульс
            signal = get_radio_pulse(varargin{:});
        case 'AM'       % тип АМ сигнал
            signal = get_AM(varargin{:});
    end
end
% Отображение сигналов
function plot_signal(varargin)
%PLOT_SIGNAL Отобразить сигнал, в случае нескольких сигналов, они
%передаются по горизонтали
% Аргументы:
%   один аргумент - значения сигнала, тогда по оси абсцис отсчеты
%   два аргумента - {временная шкала или шаг дискретизации, значения амплитуд}
%   три аргумента - третьим передается название окна
    markers = {'-r','--g',':b'};
    if or(nargin>3,nargin<1); error('Неверный формат ввода'); end
    switch nargin
        case 1
            signal = varargin{1};
            x_axis = 0:length(signal)-1;
        case {2;3}
            signal = varargin{2};
            if length(varargin{1})==1 
                x_axis = 0:varargin{1}:varargin{1}*(length(signal)-1);
            else;   x_axis = varargin{1}; 
            end
            if nargin==2; tab_name='Сигнал';
            else; tab_name=varargin{3}; end
    end
    figure('Name', tab_name,'WindowStyle','docked','Color','white')
    hold on;
    for i=1:size(signal,1)
        plot(x_axis*1e6,signal(i,:),markers{i})
    end
    xlabel('t, мкc');
end
% Отображение спектр сигнала
function [varargout] = plot_spectum(t_d,signal)
    t_d_us = t_d*1e6;
    [N_signals, n_signal] = size(signal);
    t_window = n_signal*t_d_us;
    f_step_mhz = 1/t_window;
    if ~mod(n_signal,2)
        f_axis = 0:f_step_mhz:n_signal*f_step_mhz/2;
        f_axis = [-f_axis(end-1:-1:2), f_axis];
    else
        warning('нечетное число отсчетов');
        f_axis = 0:f_step_mhz:(n_signal-1)*f_step_mhz/2;
        f_axis = [-f_axis(end:-1:2), f_axis];
    end
    signal_spectrum = 1/n_signal*fft(signal,n_signal,2);
    h = zeros(1,N_signals);
    figure('Name','Cпектр', 'WindowStyle','docked', 'Color','white')
%     hold on;
    title(['Частота дискретизации ' num2str(1/t_d_us) ' МГц'])
    for i=1:N_signals
        subplot(N_signals, 1, i)
        h(i) = plot(f_axis,fftshift(abs(signal_spectrum(i,:)))); %,markers{i});
        xlabel('Частота, МГц')
    end  
    if nargout > 0; varargout{1} = h; end
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

%% Функции формирования фильтров
% Применить к сигналу фильтр с "идеальной" АЧХ
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
% Применить к сигналу фильтр Баттерворта
function signal_out = apply_butt_filter(t_d, filter_type, f_cut_hz, filter_order, signal_in)
% Ф-я фильтрации фильтром с АЧХ Баттерворта
% Аргументы:
%   t_d - интервал дискретизации,с; 
%   filter_type: 'LP' - НЧ, 'HP' - ВЧ, 'BP' - полосовой, 'S' - заградительный
%   f_cut_hz - частота среза, Гц;
%   filter_order - порядок фильтра
%   signal_in - входной сигнал (вектор строка)
    [b, a] = create_butt_filter(t_d, filter_type, f_cut_hz, filter_order);  % формирование коэффициентов фильтра Баттерворта
    signal_out = filter(b,a,signal_in.').';     % фильтрация сигнала
end
% Применить к сигналу фильтр Чебышёва
function signal_out = apply_cheb2_filter(t_d, filter_type, f_cut_hz, filter_order, bnd_att_db, signal_in)
% Ф-я фильтрации фильтром с АЧХ Чебышева
% Аргументы:
%   t_d - интервал дискретизации,с; 
%   filter_type: 'LP' - НЧ, 'HP' - ВЧ, 'BP' - полосовой, 'S' - заградительный
%   f_cut_hz - частота среза, Гц;
%   filter_order - порядок фильтра;
%   bnd_att_db - подавление боковых лепестков АЧХ;
%   signal_in - входной сигнал (вектор строка)
    [b,a] = create_cheb2_filter(t_d, filter_type, f_cut_hz, filter_order, bnd_att_db);  % формирование фильтра
    signal_out = filter(b,a,signal_in.').';                 % фильтрация сигнала
end
% Расчет коэффициентов фильтра Баттерворта
function [b, a] = create_butt_filter(t_d, filter_type, f_cut_hz, filter_order)
    % Ф-я расчета коэффициентов передаточной функции фильтра Баттерворта
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
end
% Расчет коэффициентов Фильтра Чебышева
function [b, a] = create_cheb2_filter(t_d, filter_type, f_cut_hz, filter_order, bnd_att_db)
% Ф-я расчета коэффициентов передаточной характеристики фильтр с АЧХ Чебышева
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
end

function [b, a] = create_cheb1_filter(t_d, filter_type, f_cut_hz, filter_order, bnd_att_db)
% Ф-я расчета коэффициентов передаточной характеристики фильтр с АЧХ Чебышева
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
    [b,a] = cheby1(filter_order, bnd_att_db, w_n, f_name);  % формирование фильтра
end