clc; clear; close all;

%% Parameters
K = 5;   % Number of users
M = 10;  % Number of APs
W = 10e6; % Bandwidth (Hz)
P_max = 1; % Max transmission power per AP (W)
sigma2 = 1e-9; % Noise power (W)
epsilon = 1e-6; % Convergence tolerance
max_iter = 100; % Max iterations for Dinkelbach's algorithm

% Generate random channel gains (Rayleigh fading)
H = abs(sqrt(0.5) * (randn(M, K) + 1i * randn(M, K))); 

% Initialize power allocation randomly
P = P_max * rand(M, K);

% Circuit and backhaul power consumption (fixed)
P_circuit = 0.1 * ones(M, 2);
P_backhaul = 0.05 * ones(M, 1);

%% Dinkelbach’s Algorithm for EE Maximization
lambda = 0; % Initial lambda
iter = 0;
while iter < max_iter
    iter = iter + 1;

    % Compute SINR for each user
    SINR = zeros(K, 1);
    for k = 1:K
        signal = sum(P(:, k) .* H(:, k).^2);
        interference = sum(sum(P .* H.^2)) - signal;
        SINR(k) = signal / (sigma2 + interference);
    end

    % Compute sum-rate
    R_sum = W * sum(log2(1 + SINR));

    % Compute total power consumption
    P_tx = sum(P(:));
    P_total = P_tx + sum(P_circuit) + sum(P_backhaul);

    % Compute EE
    EE = R_sum / P_total;

    % Convergence check
    if abs(EE - lambda) < epsilon
        break;
    end

    % Update lambda
    lambda = EE;

    % Solve convex subproblem: Power allocation update (using water-filling)
    for m = 1:M
        for k = 1:K
            P(m, k) = min(P_max, (1 / log(2)) * (1 / (lambda * (sigma2 + sum(P(:, k) .* H(:, k).^2)))));
        end
    end
end

%% Display Results
disp(['Optimized EE: ', num2str(EE), ' bits/Joule']);
disp(['Total Power Consumption: ', num2str(P_total), ' W']);
disp(['Total Sum-Rate: ', num2str(R_sum / 1e6), ' Mbps']);
