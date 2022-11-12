% �������� ���������
f_d = 20e6;                          % ������� �������������, ��
t_d = 1/f_d;                        % ������ �������������
t_impulse = 20e-6;                % ������������ �������� (t chip)
t_window = 80e-6;                % ������������ ���������������� ���������
f_carrier = 1.80e6;                  % ������� �������, ��
f_mod = 0.6e6; i_mod = 0.8; % ��������� ��������� (�������, �������)
%% ������ �3 ������������� ��������
% ������ - ��������� ������ ������� � �������������� ������������ � �������-����������
% ������������ ��������� ���������
% s_c = generate_single_chip('chirp', t_d, t_window, t_impulse, f_carrier, f_mod);    % ��� �������
s_v = generate_single_chip('video', t_d, t_window, t_impulse);                      % �/�������
% s_r = generate_single_chip('radio', t_d, t_window, t_impulse, f_carrier);           % �/�������
% s_am = generate_single_chip('AM', t_d, t_window, f_carrier, f_mod, i_mod);          % �� ������
%% ����������� ��������
signal = circshift(s_v,floor(length(s_v)/2));  % ������� ������ ��� �����������
plot_signal(t_d,signal)     % ���������� ������
plot_spectum(t_d,signal);   % ���������� ������ �������
%% ��������� ��������
% �������� ���������
signal_in = signal;         % ������� ������ ��� ���������� ���������
f_type = 'HP';              % ��� ������� (LP/HP/BP/S)
f_cut_hz = 1/t_impulse;             % ������� �����, ��
but_order = 3;  cheb_order = 3;     % ������� ��������
s_band_attenuation_db = 20;         % ����������� ���������� ��������� ��� ������� ��������

signal_out_idl = apply_ideal_filter(t_d, f_type, f_cut_hz, signal_in);              % ��������� �.
signal_out_btr = apply_butt_filter(t_d, f_type, f_cut_hz, but_order, signal_in);    % �. �����������
signal_out_chb = apply_cheb2_filter(t_d, f_type, f_cut_hz, cheb_order, ...
    s_band_attenuation_db, signal_in);                                              % �. ��������
% ����������� �����������
signal_out = [signal_out_idl; signal_out_btr; signal_out_chb];  % ������������ ���-����� ����������
plot_signal(t_d,signal_out)                 % ���������� ������������� �������
legend('ideal','butterworth','chebyshev')   % ��������� �������
hndls = plot_spectum(t_d,signal_out);       % ���������� ������ (�-� ���������� ��������� �� ������)
% ��������� �������
legend(hndls(1),'ideal');
legend(hndls(2),'butterworth');
legend(hndls(3),'chebyshev');
%% ������ �4 ��������-��������� ������
% ������������ ������������������ ��������
n_chips = 8;
s_c = generate_sequence('chirp', t_d, n_chips, t_impulse, f_carrier, f_mod);   % ��� �������
s_r = generate_sequence('radio', t_d, n_chips, t_impulse, f_carrier);   % �/�������
s_am = generate_sequence('AM', t_d, n_chips, t_impulse, f_carrier, f_mod);     % �� ������
plot_signal(t_d,[s_c;s_r;s_am])

%% ����������� ������������������ ��� ���������
figure('Name', '��� ���','WindowStyle','docked','Color','white')
subplot(3,1,1); spectrogram(s_c,6,2,'yaxis');       % ������������� (�������� ����)
subplot(3,1,2); spectrogram(s_c,100,80,'yaxis');    % ������������� (������� ����, �������� ���)
subplot(3,1,3); spectrogram(s_c,100,0,'yaxis');     % ������������� (������� ����)
figure('Name', '����. ���','WindowStyle','docked','Color','white')
cwt(s_c,f_d,'amor');

%% ����������� ������������������ ��������������
figure('Name', '��� �.���.','WindowStyle','docked','Color','white')
subplot(3,1,1); spectrogram(s_r,6,2,n_fft,'yaxis');     % ������������� (�������� ����)
subplot(3,1,2); spectrogram(s_r,400,390,400,'yaxis');   % ������������� (������� ����, �������� ���)
subplot(3,1,3); spectrogram(s_r,400,0,400,'yaxis');     % ������������� (������� ����)
figure('Name', '����. �.���','WindowStyle','docked','Color','white')
cwt(s_r,f_d,'amor');

%% ����������� ������������������ �� �����
figure('Name', '��� �� ���','WindowStyle','docked','Color','white')
subplot(3,1,1); spectrogram(s_am,6,2,'yaxis');      % ������������� (�������� ����)
subplot(3,1,2); spectrogram(s_am,400,80,'yaxis');   % ������������� (������� ����, �������� ���)
subplot(3,1,3); spectrogram(s_am,400,0,'yaxis');    % ������������� (������� ����)
figure('Name', '����. �� ���','WindowStyle','docked','Color','white')
cwt(s_am,f_d);

%% ��������� ������� ������ ����������
% ������������ �������� (����������� �������)
function signal = generate_single_chip(type,varargin)
    switch type     % ��������� ������� ������������ ������� � ����������� �� ��������� ����
        case 'video'    % ��� �/�������
            signal = get_video_pulse(varargin{:});
        case 'radio'    % ��� �/�������
            signal = get_radio_pulse(varargin{:});
        case 'chirp'
            signal = get_chirp_pulse(varargin{:});
        case 'AM'       % ��� �� ������
            signal = get_AM(varargin{:});
    end
end
% ������������ ������������������ �����
function signal = generate_sequence(type,t_d, n_chips, t_imp, f_low, varargin)
%�-� ������������ ������ ������������������ �������� ��������������� ��
%�������
% ���������:
%   type - ��� �������
%   t_d - �������� �������������
%   n_chips - ����� ���������
%   t_imp - ������������ ������ ��������
%   f_low - ������ �������
% �������������� ���������:
%   z������ ���������
    switch type     % ��������� ������� ������������ ������� � ����������� �� ��������� ����
        case 'chirp'    % ��� �/�������
            f_mod =  varargin{1};
            get_chip = @(i) get_chirp_pulse(t_d, t_imp, t_imp, f_low, f_mod);
        case 'radio'    % ��� �/�������
            random_freq = f_low + randi(5,1,n_chips).*8/t_imp;       % ��������� ����������� ��� �/��������
            get_chip = @(i) get_radio_pulse(t_d, t_imp, t_imp, random_freq(i));
        case 'AM'       % ��� �� ������
            f_mod = varargin{1};   % ���������� ������� ��������� 1/10 ������� �������
            random_freq = f_low + (3*f_mod)*randi(5,1,n_chips);    % ������������ ��������� ������
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
% ����������� ��������
function [varargout] = plot_signal(varargin)
%PLOT_SIGNAL ���������� ������, � ������ ���������� ��������, ���
%���������� �� �����������
% ���������:
%   ���� �������� - �������� �������, ����� �� ��� ������ �������
%   ��� ��������� - {��������� ����� ��� ��� �������������, �������� ��������}
    switch nargin
        case 1          % ���� ������� ������ ������
            signal = varargin{1};
            x_axis = 0:length(signal)-1;    % ��������� ��� � ��������
        case 2          % ���� �������� ��������� ���������, ������ � 
            signal = varargin{2};
            if length(varargin{1})==1       % ���� ������� �������� �������������
                x_axis = 0:varargin{1}:varargin{1}*(length(signal)-1);
            else;   x_axis = varargin{1};   % ���� �������� ��� �������
            end
        otherwise
            error('�������� ������ �����: ����� ����������');
    end
    h = zeros(size(signal,1));              % ������ � ����������� �� �������
    markers = {'-r','--g',':b','-.k','-c','--m',':y','-.b'};    % ��������� �����
    figure('Name','������', 'WindowStyle','docked', 'Color','white'); hold on;
    for i=1:size(signal,1)
        h(i) = plot(x_axis*1e6,signal(i,:),markers{i});         % ������ � ������. ���������
    end
    xlabel('t, ��c');
    if nargout==1; varargout{1}=h; end  % ���� �������������� �����
end
% ����������� ������ �������
function [varargout] = plot_spectum(t_d,signal)
% ����������� �������� � ����� ���� �� �������� ���� ��� ������
% ���������:
%   t_d - �������� �������������
%   signal - ������� ������� (���������������)
    t_d_us = t_d*1e6;
    [N_signals, n_signal] = size(signal);
    n_fft = 10*n_signal;
    t_window_us = n_signal*t_d_us;
    f_step_mhz = 1/t_window_us*n_signal/n_fft;
    if ~mod(n_signal,2)
        f_axis = 0:f_step_mhz:n_fft*f_step_mhz/2;
        f_axis = [-f_axis(end-1:-1:2), f_axis];
    else
        warning('�������� ����� ��������');
        f_axis = 0:f_step_mhz:(n_fft-1)*f_step_mhz/2;
        f_axis = [-f_axis(end:-1:2), f_axis];
    end
    signal_spectrum = 1/n_signal*fft(signal,n_fft,2);
    h = zeros(1,N_signals);
    figure('Name','C�����', 'WindowStyle','docked', 'Color','white')
    title(['������� ������������� ' num2str(1/t_d_us) ' ���'])
    for i=1:N_signals
        subplot(N_signals, 1, i)
        h(i) = plot(f_axis,fftshift(abs(signal_spectrum(i,:))));
        xlabel('�������, ���')
    end  
    if nargout==1; varargout{1} = h; end
end

%% ������� ������������ ����
% ������������ �/��������
function signal = get_video_pulse(t_d, t_window, t_imp)
% �-� ������������ �/��������
% ���������:    
%   % t_d - ������ �������������, �; t_window - ������ ������������, �
%   t_imp - ������������ ��������, �
    n_signal = floor(t_window/t_d); % ����� �������� � ��������������� ������.
    n_imp = floor(t_imp/t_d);       % ����� �������� ��������
    signal = zeros(1,n_signal);     % ������������ ������ �������
    signal(1:n_imp)=1;              % ��������� ������ n_imp ��������
end
% ������������ �/��������
function signal = get_radio_pulse(t_d, t_window, t_imp, f_carrier_hz)
% �-� ������������ �/��������
% ���������:
%   t_d - �������� �������������, �;    t_window - �������� ������������, �
%   t_imp - ������������ ��������, �;   f_carrier_hz - ������� �������, ��
    n_signal = floor(t_window/t_d); % ����� �������� � ��������������� ������.
    n_imp = floor(t_imp/t_d);       % ����� �������� ��������
    pulse_amp = zeros(1, n_signal);     pulse_amp(1:n_imp)=1; % ������������ ����.
    t_axis = linspace(0, t_window, n_signal);   % ������������ ��������� ���
    carrier = sin(2*pi*f_carrier_hz*t_axis);    % ������������ �������
    signal = pulse_amp.*carrier;    % ������� �� ����������� ��������
end
% ������������ ��/�������
function signal = get_AM(t_d, t_window, f_carrire_hz, f_mod_hz, i_mod)
% �-� ������������ �� �������
%   t_d - �������� �������������, �;    t_window - �������� ������������, �
%   f_carrier_hz - ������� �������, ��;
%   f_mod_hz - ������� ���������, ��;   i_mod - ������� ��������� (Amin/Amax)
    n_signal = floor(t_window/t_d);                 % ����� �������� � ��������� �������.
    t_axis = linspace(0, t_window, n_signal);       % ��� �������
    am_mult = 1/(2/i_mod); am_shift = 1-am_mult;    % ��������� ��� ������� ��������� ��
    a_modulation = sin(2*pi*f_mod_hz*t_axis)*am_mult+am_shift;  % ���������
    signal = sin(2*pi*f_carrire_hz*t_axis) .* a_modulation;     % ������������ �� �������
end
% ������������ ���
function signal = get_chirp_pulse(t_d, t_window, t_imp, f_start_hz, f_chirp_hz)
% �-� ������������ �/��������
% ���������:
%   t_d - �������� �������������, �;    t_window - �������� ������������, �
%   t_imp - ������������ ��������, �;   f_center_hz - ���������� �������, ��
%   f_chirp_hz - ������ ������, ��
    n_signal = floor(t_window/t_d);     % ����� �������� � ��������������� ������.
    n_imp = floor(t_imp/t_d);           % ����� �������� ��������
    t_axis = linspace(0, t_imp, n_imp);   % ������������ ��������� ���
    chirp_band = f_chirp_hz/t_imp;
    chirp_sin_arg = f_start_hz + chirp_band/2*t_axis;
    signal = zeros(1,n_signal);
    signal(1:n_imp) = sin(2*pi*chirp_sin_arg.*t_axis);     % ������������ �������
%     signal(n_imp:end) =  sin(2*pi*f_moment(n_imp).*t_axis(n_imp:end));     
    
end

%% ������� ������������ ��������
% ������������ ������� � ��������� ��� � ��������(-���) ��������� �����
function signal_out = apply_ideal_filter(t_d, filter_type, f_cut_hz, signal_in)
% �-� ���������� �������� � ��������� ���
% ���������:
%   t_d - �������� �������������,�; 
%   filter_type: 'LP' - ��, 'HP' - ��, 'BP' - ���������, 'S' - ��������������
%   f_cut_hz - ������� �����, ��;
%   signal_in - ������� ������ (������ ������)
    f_d = 1/t_d;    % ������ ������� �������������
    n_signal = length(signal_in);
    f_axis = linspace(0,f_d/2,n_signal/2);
    if length(f_cut_hz)==1         % ���� ���� �������� ������� �����, �� �� ��� �� ������
        [~, f_cut_i] = find(f_axis>=f_cut_hz,1);        % ����������� ������ �������
        filter_f_char = ones(1,ceil((n_signal+1)/2));   % ��� 2n -> n+1; ��� 2n+1 -> n+1
        if      strcmp(filter_type,'HP'); filter_f_char(1:f_cut_i-1) = 0;   % ���� �� 
        elseif  strcmp(filter_type,'LP'); filter_f_char(f_cut_i:end) = 0;   % ���� ��
        else;   error('������ ������: ��� �������'); 
        end
        % ������������ ��� ��� �������������� �������
        filter_f_char = [...
            filter_f_char, ...
            filter_f_char(end-(mod(n_signal,2)==0):-1:2)]; % � �������������� �������
    elseif length(f_cut_hz)==2     % ���� ��� �������� ������� �����, �� �� ��� �� ������
        %TODO: ����������� �������� �� � ��
    else; error('������ �����: ��������� ������� �������')
    end
    signal_in_sp = fft(signal_in, n_signal, 2);     % ��� �������
    signal_out_sp = signal_in_sp .* filter_f_char;  % ����������
    signal_out = ifft(signal_out_sp, n_signal, 2);  % �������� ��� ������� �������
    % �������� ��������� �������: ���� �� �����������, �� �������� ������
    if ~isreal(signal_out)  
        warning('�������� ������ ��� ������������ �������: ������ ����� ���������� �����������')
    end
end
% ������������ ������� ����������� � ������� �������� ����� � ��������
function signal_out = apply_butt_filter(t_d, filter_type, f_cut_hz, filter_order, signal_in)
% �-� ���������� �������� � ��� �����������
% ���������:
%   t_d - �������� �������������,�; 
%   filter_type: 'LP' - ��, 'HP' - ��, 'BP' - ���������, 'S' - ��������������
%   f_cut_hz - ������� �����, ��;
%   filter_order - ������� �������
%   signal_in - ������� ������ (������ ������)
    if length(f_cut_hz)==1      % ���� ���� �������� - �� ��� �� ������
        switch filter_type
            case 'LP'; f_name = 'low';  % ������ ��
            case 'HP'; f_name = 'high'; % ������ ��
            otherwise; error('������ ������: ��� �������'); 
        end
    elseif length(f_cut_hz)==2  % ���� ��� �������� ������� - �� �� ��
        switch filter_type
            case 'BP'; f_name = 'bandpass';     % ��������� ������
            case 'S'; f_name = 'stop';          % �������������� ������
            otherwise; error('������ ������: ��� �������'); 
        end
    else; error('������ �����: ��������� ������� �������');
    end
    f_d = 1/t_d;            % ������� �������������
    w_n = f_cut_hz/(f_d/2); % ������������� ������� �����
    [b,a] = butter(filter_order, w_n, f_name);  % ������������ ������������� ������� �����������
    signal_out = filter(b,a,signal_in.').';     % ���������� �������
end
% ������������ ������� �������� � �������� �������� ����� � ��������
% (����������� ���������)
function signal_out = apply_cheb2_filter(t_d, filter_type, f_cut_hz, filter_order, bnd_att_db, signal_in)
% �-� ���������� �������� � ��� �����������
% ���������:
%   t_d - �������� �������������,�; 
%   filter_type: 'LP' - ��, 'HP' - ��, 'BP' - ���������, 'S' - ��������������
%   f_cut_hz - ������� �����, ��;
%   filter_order - ������� �������;
%   bnd_att_db - ���������� ������� ��������� ���;
%   signal_in - ������� ������ (������ ������)
    if length(f_cut_hz)==1  % ���� �������� ���� ������� ������� �����
        switch filter_type
            case 'LP'; f_name = 'low';  % �� ������
            case 'HP'; f_name = 'high'; % �� ������
            otherwise; error('������ ������: ��� �������'); 
        end
    elseif length(f_cut_hz)==2  % ���� �������� ��� �������� ������ �����
        switch filter_type
            case 'BP'; f_name = 'bandpass'; % ��������� ������
            case 'S'; f_name = 'stop';      % �������������� ������
            otherwise; error('������ ������: ��� �������'); 
        end
    else; error('������ �����: ��������� ������� �������');
    end
    f_d = 1/t_d;            % ������� �������������, ��
    w_n = f_cut_hz/(f_d/2); % ������������� ������� ����� (w/w_d)
    [b,a] = cheby2(filter_order, bnd_att_db, w_n, f_name);  % ������������ �������
    signal_out = filter(b,a,signal_in.').';                 % ���������� �������
end
