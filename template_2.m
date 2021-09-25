
time_window_us = 90;    % âðåìåííîé èíòåðâàë íà êîòîðîì ðàññìàòðèâàåòñÿ êîððåëÿöèÿ, ìêñ
sigma_normal = 24;      % ÑÊÎ ãàóññîâñêîãî ñëó÷àéíîãî ïðîöåññà
corr_int_us = 10;        % èíòåðàâàë êîððåëÿöèè, ìêñ
n_samples = 32;         % ÷èñëî óñðåäíÿåìûõ ðåàëèçàöèé
f_sampling_mhz = 20;    % ÷àñòîòà äèñêðåòèçàöèè, ÌÃö
n_counts = ceil(time_window_us*f_sampling_mhz); % ÷èñëî îòñ÷åòîâ

%% Ôîðìèðîâàíèå êîððåëèðîâàííîãî ñèãíàëà
x = sigma_normal * randn(n_samples, n_counts);          % ðåàëèçàöèÿ íîðìàëüíîãî íåêîððåëèðîâàííîãî øóìà
y = corr_filter(x, time_window_us, corr_int_us);          % ðåàëèçàöèè êîððåëèðîâàííîãî íîðì. øóìà
rout = zeros(n_samples, 2*n_counts-1);
for i=1:n_samples
    rout(i,:) = xcorr(y(i,:), 'coeff'); % ðàñ÷åò êîððåëÿöèè äëÿ êàæäîé ðåàëèçàöèè
end
rout = mean(rout,1);        % óñðåäíåíèå ðåàëèçàöèé

%% Ãðàôèêà
t_axis_us = linspace(0,time_window_us,n_counts);              % îñü âðåìåíè äëÿ ðåàëèçàöèé
shift_axis_us = ...
    linspace(-time_window_us, time_window_us, 2*n_counts-1);    % îñü âðåìåíè äëÿ êîððåëÿöèîííûõ ôóíêöèé
anal_corr_fun = exp(-abs(shift_axis_us)/(corr_int_us/2));   % àíàëèò. êîððåëÿöèîííàÿ ôóíêöèÿ äëÿ ãðàôèêà

% ãðàôèêè äëÿ îäíîé ðåàëèçàöèè
figure('Name', 'single sample', 'color', 'white','WindowStyle','docked');
subplot(2,1,1);
plot(t_axis_us, y(1,:));    % âûâîä ðåàëèçàöèé ðåçóëüòèðóþùåãî êîððåëèðîâàííîãî ïðîöåññà
subplot(2,1,2);
plot(shift_axis_us, xcorr(y(1,:),'coeff'));  % 

% ãðàôèêè óñðåäíåííûõ êîððåëÿöèîííûõ ôóíêöèé
figure('Name', 'averaged', 'color', 'white','WindowStyle','docked');
grid on; hold on; 
plot(shift_axis_us,rout);
plot(shift_axis_us,anal_corr_fun);
legend('óñðåäíåííàÿ ÊÔ','')
xlabel('\tau'); ylabel('R(\tau');

%% Ôèëüòð, îáåñïå÷èâàþùèé êîððåëÿöèþ îòñ÷åòîâ
function [f] = corr_filter(Z,rL_us,corr_us)
    [~, signal_length] = size(Z);                   % îïåðåäåëåíèå 
    x = linspace(-rL_us/2,rL_us/2,signal_length);   % ïåðåñ÷åò îòñ÷åòîâ âî âðåìÿ
    F = exp(-abs(x)/(corr_us/2));                   % êîððåëÿöèîííàÿ ôóíêöèÿ
    f = sqrt(2)*ifft(fft(Z,signal_length,2).*sqrt(fft(F,signal_length,2)),signal_length,2); % ðåçóëüòàò ïðåîáðàçîâàíèÿ
end

