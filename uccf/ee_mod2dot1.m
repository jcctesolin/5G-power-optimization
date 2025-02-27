clc; clear; close all;

%% Parameters
K = 2 ;   % Number of users
L = 4 ;  % Number of APs
W = 10e6; % Bandwidth (Hz)
P_max = 5; % Max transmission power per AP (W)
sigma2 = 1e-9; % Noise power (W)
epsilon = 1e-6; % Convergence tolerance
max_iter = 100; % Max iterations for Dinkelbach's algorithm
nsrv = 4;
ngh_size = 6;

%% Generate random channel gains (Rayleigh fading)
%H = abs(sqrt(0.5) * (randn(L, K) + 1i * randn(L, K))); 
H = [0.5 0.8 1.0 0.6; 
    1.0 0.5 0.4 0.6]';

[sorted_vals, sorted_indices] = sort(H);
%[row, col] = ind2sub(size(H), sorted_indices);

%% Creating clusters

%D = zeros(L);
clusterset = cell(1:K);

for k=1:K
    cluster = [0 0 0 0];
    cluster_base = [0 0 0 0];
    for j=1:nsrv
        if j==1
            cluster(1,sorted_indices(j,k))=1;
            cluster_base = cluster;
        else
            cluster = cluster_base;
            cluster(1,sorted_indices(j,k))=1;
        end
        clusterset{k}{end+1}=cluster;
    end    
end

%% Initialize power allocation randomly
%P = P_max * rand(L, K);
P = P_max * ones(L, K);
%P = D * P;
%P = P_max * D;

% Circuit and backhaul power consumption (fixed)
%P_circuit = 0.1 * ones(L, 1);
P_circuit = 0;
%P_backhaul = 0.05 * ones(L, 1);
P_backhaul = 0;

%% Select best cluster set among UEs
bag = {};
for k=1:K
    bag{end+1} = clusterset{K}{randi(length(clusterset{K}))};
end

%% Dinkelbach’s Algorithm for EE Maximization
lambda = 0; % Initial lambda
iter = 0;
while iter < max_iter
    iter = iter + 1;

    % Compute SINR for each user
    SINR = zeros(K, 1); % zero matrix, matrix initialization
    for k = 1:K
        %signal = sum(P(:, k) .* H(:, k).^2);
        %signal = sum(P(:, k) .* D.* H(:, k).^2);
        signal =  bag{k}*(H(:, k).^2.* P(:, k)) ;
        %interference = sum(sum(P .* H.^2)) - signal;
        interference = sum(H(:, k).^2.* P(:, k)) - signal ;
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
    for l = 1:L
    %for l = cluster
        for k = 1:K
            %P(l, k) = min(P_max, (1 / log(2)) * (1 / (lambda * (sigma2 + sum(P(:, k) .* H(:, k).^2)))));
            P(l, k) = min(P_max, (1 / log(2)) * (1 / (lambda * (sigma2 + sum(P(:, k) .* H(:, k).^2)))));
        end
    end
end

%% Display Results
disp(['Optimized EE: ', num2str(EE), ' bits/Joule']);
disp(['Total Power Consumption: ', num2str(P_total), ' W']);
disp(['Total Sum-Rate: ', num2str(R_sum / 1e6), ' Mbps']);
disp(['Num Iter: ', num2str(iter),' iterations']);
