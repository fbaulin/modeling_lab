
time_window_us = 90;    % временной интервал на котором рассматривается корреляция, мкс
sigma_normal = 24;      % СКО гауссовского случайного процесса
corr_int_us = 10;        % интеравал корреляции, мкс
n_samples = 32;         % число усредняемых реализаций
f_sampling_mhz = 20;    % частота дискретизации, МГц
n_counts = ceil(time_window_us*f_sampling_mhz); % число отсчетов

%% Формирование коррелированного сигнала
x = sigma_normal * randn(n_samples, n_counts);          % реализация нормального некоррелированного шума
y = corr_filter(x, time_window_us, corr_int_us);          % реализации коррелированного норм. шума
rout = zeros(n_samples, 2*n_counts-1);
for i=1:n_samples
    rout(i,:) = xcorr(y(i,:), 'coeff'); % расчет корреляции для каждой реализации
end
rout = mean(rout,1);        % усреднение реализаций

%% Графика
t_axis_us = linspace(0,time_window_us,n_counts);              % ось времени для реализаций
shift_axis_us = ...
    linspace(-time_window_us, time_window_us, 2*n_counts-1);    % ось времени для корреляционных функций
anal_corr_fun = exp(-abs(shift_axis_us)/(corr_int_us/2));   % аналит. корреляционная функция для графика

% графики для одной реализации
figure('Name', 'single sample', 'color', 'white','WindowStyle','docked');
subplot(2,1,1);
plot(t_axis_us, y(1,:));    % вывод реализаций результирующего коррелированного процесса
subplot(2,1,2);
plot(shift_axis_us, xcorr(y(1,:),'coeff'));  % 

% графики усредненных корреляционных функций
figure('Name', 'averaged', 'color', 'white','WindowStyle','docked');
grid on; hold on; 
plot(shift_axis_us,rout);
plot(shift_axis_us,anal_corr_fun);
legend('усредненная КФ','')
xlabel('\tau'); ylabel('R(\tau');

%% Фильтр, обеспечивающий корреляцию отсчетов
function [f] = corr_filter(Z,rL_us,corr_us)
    [~, signal_length] = size(Z);                   % опеределение 
    x = linspace(-rL_us/2,rL_us/2,signal_length);   % пересчет отсчетов во время
    F = exp(-abs(x)/(corr_us/2));                   % корреляционная функция
    f = sqrt(2)*ifft(fft(Z,signal_length,2).*sqrt(fft(F,signal_length,2)),signal_length,2); % результат преобразования
end

