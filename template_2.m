time_window_us = 90;    % временной интервал на котором рассматриваетс€ коррел€ци€, мкс
sigma_normal = 24;      % — ќ гауссовского случайного процесса
corr_int_us = 10;        % интеравал коррел€ции, мкс
n_samples = 32;         % число усредн€емых реализаций
f_sampling_mhz = 20;    % частота дискретизации, ћ√ц
n_counts = ceil(time_window_us*f_sampling_mhz); % число отсчетов

%% ‘ормирование коррелированного сигнала
x = sigma_normal * randn(n_samples, n_counts);          % реализаци€ нормального некоррелированного шума
y = corr_filter(x, time_window_us, corr_int_us);          % реализации коррелированного норм. шума
rout = zeros(n_samples, 2*n_counts-1);
for i=1:n_samples
    rout(i,:) = xcorr(y(i,:), 'coeff'); % расчет коррел€ции дл€ каждой реализации
end
rout = mean(rout,1);        % усреднение реализаций

%% √рафика
t_axis_us = linspace(0,time_window_us,n_counts);              % ось времени дл€ реализаций
shift_axis_us = ...
    linspace(-time_window_us, time_window_us, 2*n_counts-1);    % ось времени дл€ коррел€ционных функций
anal_corr_fun = exp(-abs(shift_axis_us)/(corr_int_us/2));   % аналит. коррел€ционна€ функци€ дл€ графика

% графики дл€ одной реализации
figure('Name', 'single sample', 'color', 'white','WindowStyle','docked');
subplot(2,1,1);
plot(t_axis_us, y(1,:));    % вывод реализаций результирующего коррелированного процесса
subplot(2,1,2);
plot(shift_axis_us, xcorr(y(1,:),'coeff'));  % 

% графики усредненных коррел€ционных функций
figure('Name', 'averaged', 'color', 'white','WindowStyle','docked');
grid on; hold on; 
plot(shift_axis_us,rout);
plot(shift_axis_us,anal_corr_fun);
legend('усредненна€  ‘','')
xlabel('Tau'); ylabel('R(Tau');

%% ‘ильтр, обеспечивающий коррел€цию отсчетов
function [f] = corr_filter(Z,rL_us,corr_us)
    [~, signal_length] = size(Z);                   % опеределение 
    x = linspace(-rL_us/2,rL_us/2,signal_length);   % пересчет отсчетов во врем€
    F = exp(-abs(x)/(corr_us/2));                   % коррел€ционна€ функци€
    f = sqrt(2)*ifft(fft(Z,signal_length,2).*sqrt(fft(F,signal_length,2)),signal_length,2); % результат преобразовани€
end

