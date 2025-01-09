function [ber_zf,ber_mmse] = MU_MIMO_DL_precoders(snr_db)

%  MU-MIMO wireless communication system operating with 
%  spatial multiplexing with linear ZF and MMSE precoders. 
%  The model describes the downlink of a system in which each 
%  user is equipped with N_U antennas and the access point 
%  has N_A antennas. The network has K users.
%
%  July 2012
%  Authors: Rodrigo C. de Lamare 
%-------------------------------------------------------------------
% Configuration 

% QPSK


% Number of users
K = 3;

% Number of antennas per user
N_U = 2;

% Number of receive antennas
N_A = 8;

% Packet size
N = 100;

% Signal-to-noise ratio (SNR) manipulation
sigma_s = 1; % signal square root power
snr = exp(snr_db*log(10)/10); % conversion from dB  
sigma_n = sqrt( (K*N_U)/snr )*sigma_s; % noise standard deviation

% ZF and MMSE precoders
W = zeros(N_A,N_U*K); % linear filter for linear detector and K users
W_sic = zeros(N_A,N_U*K); % linear filter for SIC detector and K users

% Counters for errors
ne_zf = 0;
ne_mmse = 0;

% Channel Matrix
sigma_c = 1; % channel power standard deviation
H = zeros(N_U*K,N_A);
H = (sigma_c*randn(K*N_U,N_A)+sqrt(-1)*(sigma_c*randn(K*N_U,N_A)))/sqrt(2); 


% Precoders
P_mmse = H'*inv(H*H' + (sigma_n^2/sigma_s^2)*eye(K*N_U)); %(K*N_U/N_A)
P_zf = H'*pinv(H*H');    
% Normalization to ensure that the precoded signal energy is the same as
% before precoding- check Joham et al., TSP, 2005.
beta_zf = sqrt(N_A/trace(pinv(H*H')));
beta_mmse = sqrt(N_A/trace(P_mmse*P_mmse'));
        
for l = 1:N
    % QPSK symbols
    s_qpsk = [sign(randint(K*N_U,1)-0.5) + sqrt(-1)*sign(randint(K*N_U,1)-0.5)];
    % Scaling
    s = s_qpsk/sqrt(2);
    % Noise    
    n = (sigma_n/sqrt(2))*(randn(K*N_U,1)+sqrt(-1)*randn(K*N_U,1));
    % Received vector with normalization
    r_zf = H*P_zf*s + (n/beta_zf);
    r_mmse = H*P_mmse*s + (n/beta_mmse);
    % Slicer for QPSK
    s_zf = sign(real(r_zf)) + sqrt(-1)*sign(imag(r_zf));
    s_mmse = sign(real(r_mmse)) + sqrt(-1)*sign(imag(r_mmse));
    % Error counting     
    % ZF 
    [ex,ey] = find(real(s_zf)~= real(s_qpsk));
    ne_zf = ne_zf + length(ex);
    [ex,ey] = find(imag(s_zf)~= imag(s_qpsk));
    ne_zf = ne_zf + length(ex);
    % MMSE
    [ex,ey] = find(real(s_mmse)~= real(s_qpsk));
    ne_mmse = ne_mmse + length(ex);
    [ex,ey] = find(imag(s_mmse)~= imag(s_qpsk));
    ne_mmse = ne_mmse + length(ex);
 end
ber_zf = ne_zf/(2*N*K*N_U);
ber_mmse = ne_mmse/(2*N*K*N_U);