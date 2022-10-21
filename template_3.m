%% ������������ ������ �3
% ������ - ������������ ������� � ��������� �� ����������
% ��������� ��� ������������ ���������� � ��������������� ��������
% ������������ �������
f_d = 1e9;              % ������� �������������, ��
t_d = 1/f_d;            % ������ �������������
t_win = 2e-7;           % ������������ ���������������� ���������
t_ch = 8e-9;            % ������������ �������� (t chirp)
f_car = 140e6;          % ������� �������, ��
f_mod = 5e6; i_mod = 0.5;   % ��������� ��������� (�������,
% ������������ ��������
s_v = generate_single_chip('video', t_d, t_win, t_ch);                 % �/�������
s_r = generate_single_chip('radio', t_d, t_win, t_ch, f_car);          % �/�������
s_am = generate_single_chip('AM', t_d, t_win, f_car, f_mod, i_mod);    % �� ������
% ����������� ������
signal = circshift(s_v,floor(length(s_r)/2));   % ������� ������ ��� �����������
plot_signal(t_d,signal)     % ���������� ������
plot_spectum(t_d,signal);   % ���������� ������ �������
%% ��������� ��������
% ��������� ��� ������
signal_in = signal;     % ������� ������ ��� ���������� ���������
f_type = 'LP';          % ��� ������� (LP/HP/BP/S)
f_cut_hz = 1/t_ch;      % ������� �����, ��
but_order = 3;  cheb_order = 5; % ������� ��������
s_band_attenuation_db = 20;     % ����������� ���������� ��������� ��� ������� ��������

signal_out_idl = apply_ideal_filter(t_d, f_type, f_cut_hz, signal_in);              % ��������� �.
signal_out_btr = apply_butt_filter(t_d, f_type, f_cut_hz, but_order, signal_in);    % �. �����������
signal_out_chb = apply_cheb2_filter(t_d, f_type, f_cut_hz, ...
    cheb_order, s_band_attenuation_db, signal_in);              % �. ��������
% ����������� ������
signal_out = [signal_out_idl; signal_out_btr; signal_out_chb];  % ������������ ���-����� ����������
plot_signal(t_d,signal_out)                 % ���������� ������������� �������
legend('ideal','butterworth','chebyshev')   % ��������� �������
hndls = plot_spectum(t_d,signal_out);       % ���������� ������ (�-� ���������� ��������� �� ������)
legend(hndls(1),'ideal');  legend(hndls(2),'butterworth');  legend(hndls(3),'chebyshev');   % ��������� �������
%% ����� ��� � ���
[b, a] = create_cheb2_filter(t_d, f_type, f_cut_hz, cheb_order, s_band_attenuation_db);      % ������������ �������
n_imp_resp = 2^7;   % ������������ ��������� (� ��������) ��� �������� �������� ���. ��������������
n_freq_resp = 2^10;  % ����� �������� ��� � ���
figure('Name', '��','WindowStyle','docked','Color','white');   
impz(b, a, n_imp_resp, f_d)     % ���������� �� (����� �������� ������� - ��. �������)
figure('Name', '���,���','WindowStyle','docked','Color','white');   
freqz(b, a, n_freq_resp, f_d)   % ���������� ��� � ��� (����� �������� ������� - ��. �������)

%% ��������� ������� ������ ����������
% ������������ �������� (����������� �������)
function signal = generate_single_chip(type,varargin)
    switch type     % ��������� ������� ������������ ������� � ����������� �� ��������� ����
        case 'video'    % ��� �/�������
            signal = get_video_pulse(varargin{:});
        case 'radio'    % ��� �/�������
            signal = get_radio_pulse(varargin{:});
        case 'AM'       % ��� �� ������
            signal = get_AM(varargin{:});
    end
end
% ����������� ��������
function plot_signal(varargin)
%PLOT_SIGNAL ���������� ������, � ������ ���������� ��������, ���
%���������� �� �����������
% ���������:
%   ���� �������� - �������� �������, ����� �� ��� ������ �������
%   ��� ��������� - {��������� ����� ��� ��� �������������, �������� ��������}
%   ��� ��������� - ������� ���������� �������� ����
    markers = {'-r','--g',':b'};
    if or(nargin>3,nargin<1); error('�������� ������ �����'); end
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
            if nargin==2; tab_name='������';
            else; tab_name=varargin{3}; end
    end
    figure('Name', tab_name,'WindowStyle','docked','Color','white')
    hold on;
    for i=1:size(signal,1)
        plot(x_axis*1e6,signal(i,:),markers{i})
    end
    xlabel('t, ��c');
end
% ����������� ������ �������
function [varargout] = plot_spectum(t_d,signal)
    t_d_us = t_d*1e6;
    [N_signals, n_signal] = size(signal);
    t_window = n_signal*t_d_us;
    f_step_mhz = 1/t_window;
    if ~mod(n_signal,2)
        f_axis = 0:f_step_mhz:n_signal*f_step_mhz/2;
        f_axis = [-f_axis(end-1:-1:2), f_axis];
    else
        warning('�������� ����� ��������');
        f_axis = 0:f_step_mhz:(n_signal-1)*f_step_mhz/2;
        f_axis = [-f_axis(end:-1:2), f_axis];
    end
    signal_spectrum = 1/n_signal*fft(signal,n_signal,2);
    h = zeros(1,N_signals);
    figure('Name','C�����', 'WindowStyle','docked', 'Color','white')
%     hold on;
    title(['������� ������������� ' num2str(1/t_d_us) ' ���'])
    for i=1:N_signals
        subplot(N_signals, 1, i)
        h(i) = plot(f_axis,fftshift(abs(signal_spectrum(i,:)))); %,markers{i});
        xlabel('�������, ���')
    end  
    if nargout > 0; varargout{1} = h; end
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

%% ������� ������������ ��������
% ��������� � ������� ������ � "���������" ���
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
% ��������� � ������� ������ �����������
function signal_out = apply_butt_filter(t_d, filter_type, f_cut_hz, filter_order, signal_in)
% �-� ���������� �������� � ��� �����������
% ���������:
%   t_d - �������� �������������,�; 
%   filter_type: 'LP' - ��, 'HP' - ��, 'BP' - ���������, 'S' - ��������������
%   f_cut_hz - ������� �����, ��;
%   filter_order - ������� �������
%   signal_in - ������� ������ (������ ������)
    [b, a] = create_butt_filter(t_d, filter_type, f_cut_hz, filter_order);  % ������������ ������������� ������� �����������
    signal_out = filter(b,a,signal_in.').';     % ���������� �������
end
% ��������� � ������� ������ ��������
function signal_out = apply_cheb2_filter(t_d, filter_type, f_cut_hz, filter_order, bnd_att_db, signal_in)
% �-� ���������� �������� � ��� ��������
% ���������:
%   t_d - �������� �������������,�; 
%   filter_type: 'LP' - ��, 'HP' - ��, 'BP' - ���������, 'S' - ��������������
%   f_cut_hz - ������� �����, ��;
%   filter_order - ������� �������;
%   bnd_att_db - ���������� ������� ��������� ���;
%   signal_in - ������� ������ (������ ������)
    [b,a] = create_cheb2_filter(t_d, filter_type, f_cut_hz, filter_order, bnd_att_db);  % ������������ �������
    signal_out = filter(b,a,signal_in.').';                 % ���������� �������
end
% ������ ������������� ������� �����������
function [b, a] = create_butt_filter(t_d, filter_type, f_cut_hz, filter_order)
    % �-� ������� ������������� ������������ ������� ������� �����������
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
end
% ������ ������������� ������� ��������
function [b, a] = create_cheb2_filter(t_d, filter_type, f_cut_hz, filter_order, bnd_att_db)
% �-� ������� ������������� ������������ �������������� ������ � ��� ��������
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
end

function [b, a] = create_cheb1_filter(t_d, filter_type, f_cut_hz, filter_order, bnd_att_db)
% �-� ������� ������������� ������������ �������������� ������ � ��� ��������
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
    [b,a] = cheby1(filter_order, bnd_att_db, w_n, f_name);  % ������������ �������
end