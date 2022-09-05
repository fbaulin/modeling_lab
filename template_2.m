time_window_us = 400;    % ��������� �������� �� ������� ����������� �������������, ���
f_sampling_mhz = 20;    % ������� �������������, ���
sigma_normal = 24;      % ��� ������������ ���������� ��������
corr_int_us = 10;       % ��������� ����������, ���
n_samples = 32;         % ����� ����������� ����������
n_counts = ceil(time_window_us*f_sampling_mhz); % ����� ��������

%% ������������ ���������������� �������
x = sigma_normal * randn(n_samples, n_counts);          % ���������� ����������� ������������������ ����
y = corr_filter(x, time_window_us, corr_int_us);        % ���������� ���������������� ����. ����
% ������ ����������
r_out_single = xcorr(y(1,:),'coeff');   % ���� ���������� ��
rout = zeros(n_samples, 2*n_counts-1);  % ������� ��� ������ ���������� ��
for i=1:n_samples
    rout(i,:) = xcorr(y(i,:), 'coeff'); % ������ ���������� ��� ������ ����������
end
rout = mean(rout,1);                    % ���������� ���������� ��

%% �������
t_axis_us = linspace(0,time_window_us,n_counts);                % ��� ������� ��� ����������
shift_axis_us = ...
    linspace(-time_window_us, time_window_us, 2*n_counts-1);    % ��� ������� ��� �������������� �������
anal_corr_fun = exp(-abs(shift_axis_us)/(corr_int_us/2));       % ������. �������������� ������� ��� �������

% ������� ��� ����� ����������
figure('Name', 'single sample', 'color', 'white','WindowStyle','docked');
subplot(2,1,1); plot(t_axis_us, y(1,:));            % ����� ���������� ��������������� ���������������� ��������
subplot(2,1,2); plot(shift_axis_us, r_out_single);  % ����� ������ �������������� �������, ��������� �� ����� ����������

% ������� ����������� �������������� �������
figure('Name', 'averaged', 'color', 'white','WindowStyle','docked');
grid on; hold on; 
plot(shift_axis_us,rout);
plot(shift_axis_us,anal_corr_fun);
legend('��.��','�')
xlabel('\tau'); ylabel('R(\tau)');

%% ������, �������������� ���������� ��������
function [f] = corr_filter(Z,rL_us,corr_us)
    [~, signal_length] = size(Z);                   % ����������� 
    x = linspace(-rL_us/2,rL_us/2,signal_length);   % �������� �������� �� �����
    f = exp(-abs(x)/(corr_us/2));                   % ������� �������������� �������
    F = sqrt(fft(f,signal_length,2));               % ��� �������
    F = F/max(abs(F));  % ���������� ��� �� ������������ ��������: ��������� �� ������ �����������
    f = real(sqrt(2)*ifft(fft(Z,signal_length,2).*F,signal_length,2)); % ������������ ������� � ���
end

