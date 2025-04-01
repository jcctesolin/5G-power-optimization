
clc; clear; close all;

% Parameters
N = 3; % Number of users
Pmax = [10; 15; 20]; % Maximum power limits (Watts)
H = rand(1, N) * 10; % Channel gains
SINR_target = [2; 3; 4]; % Target SINR values
noise = 1e-9;

% Optimization problem
cvx_begin quiet
    variable P(N)
    minimize(sum(P)) % Objective: Minimize total power
    subject to
        P >= 0
        P <= Pmax
        H * P / (noise + sum(H * P) - H * P) >= SINR_target
cvx_end

% Display results
disp('Optimized Power Allocation:');
disp(P);


