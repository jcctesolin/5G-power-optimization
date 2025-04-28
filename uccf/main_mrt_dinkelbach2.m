% main_mrt_dinkelbach.m
% Algoritmo de Dinkelbach para maximizar eficiência energética com MRT

clear; clc;

load('channel_data_mimo.mat');  % Carregar canais

K = size(h_channels_real, 1);
L = size(h_channels_real, 2);
M = size(h_channels_real, 3);

% Recriar canais complexos
h = complex(h_channels_real, h_channels_imag);

% Parâmetros
noise_power = 1e-9;
Pc = 1;                      % Potência de circuito
Pmax = 1;                    % Potência máxima por AP
Pmin = 0.01;                 % Potência mínima por AP
Rmin = 100 / log2(exp(1));   % Taxa mínima em nats/s

% Inicialização
lambda = 0;
tolerance = 1e-3;
max_iter = 50;
iter = 0;

% Inicialização de potências
p = ones(K, L);  % potência inicial para cada usuário k em cada AP l

while iter < max_iter
    iter = iter + 1;

    % Função objetivo para dado lambda: maximize sum_rate - lambda * total_power
    fun = @(p_vec) - (sum_rate_mrt(reshape(p_vec, K, L), h, noise_power) - ...
                      lambda * (total_power_mrt(reshape(p_vec, K, L)) + Pc));

    % Restrições
    lb = Pmin * ones(K*L, 1);
    ub = Pmax * ones(K*L, 1);

    nonlcon = @(p_vec) rate_constraints(reshape(p_vec, K, L), h, noise_power, Rmin);

    %options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
    options = optimoptions('fmincon','Display','final','Algorithm','sqp');
    [p_opt_vec, ~] = fmincon(fun, p(:), [], [], [], [], lb, ub, nonlcon, options);

    p_opt = reshape(p_opt_vec, K, L);
    num = sum_rate_mrt(p_opt, h, noise_power);
    denom = total_power_mrt(p_opt) + Pc;
    new_lambda = num / denom;

    if abs(new_lambda - lambda) < tolerance
        break;
    end
    lambda = new_lambda;
end

%fprintf('Eficiência energética ótima: %.4f bps/W', lambda);
%fprintf('SumRate: %.4f bps', num);
%fprintf('TotalPower: %.4f W', denom);

disp(['Optimized EE: ', num2str(lambda/1e6,'%.2e'), ' Mbits/Joule']);
disp(['Total Power Consumption: ', num2str(denom), ' W']);
disp(['Total Sum-Rate: ', num2str(num), ' Mbps']);
disp('Power Matrix:');
disp(p_opt);