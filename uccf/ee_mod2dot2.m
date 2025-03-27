clc; clear; close all;

%% Parameters
K = 2 ;   % Number of users
L = 4 ;  % Number of APs
W = 10e6; % Bandwidth (Hz)
P_max = 10; % Max transmission power per AP (W)
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
%P = P_max * ones(L, K) * 1/K;
%P = D * P;
%P = P_max * D;

% Circuit and backhaul power consumption (fixed)
%P_circuit = 0.1 * ones(L, 1);
P_circuit = 0;
%P_backhaul = 0.05 * ones(L, 1);
P_backhaul = 0;

%% Select best cluster set among UEs

T = combinations(clusterset{:});
bag = T{:,:};
%disp(bag{1,1})

% for k=1:K
%     bag{end+1} = clusterset{K}{randi(length(clusterset{K}))};
% end

%% Dinkelbach’s Algorithm for EE Maximization

bagrows = size(bag, 1);

for bagiter=1:bagrows
    lambda = 0; % Initial lambda
    iter = 0;
    P = P_max * ones(L, K) * 1/K;
    while iter < max_iter
        iter = iter + 1;
        % Compute SINR for each user
        SINR = zeros(K, 1);
        for k = 1:K
            %signal =  bag{bagiter,k}*(H(:, k).^2.* P(:, k)) ;
            signal =  sum(bag{bagiter,k}*(H(:, k).^2.* P(:, k))) ;
            %interference = sum(sum(P .* H.^2)) - signal;
            interference = sum(H(:, k).^2.* P(:, k)) - signal ;
            SINR(k) = signal / (sigma2 + interference);
            %R(k) = W * log2(1 + SINR(k)) ;
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
        
        %dummy= sum(bag{bagiter,k}*(H(:, k).^2.* P(:, k)))/(sigma2+(sum(H(:, k).^2.* P(:, k)) - sum(bag{bagiter,k}*(H(:, k).^2.* P(:, k)))));
        %dummy_num=sum(bag{bagiter,k}*(H(:, k).^2.* P(:, k)));
        %dummy_den=(sigma2+(sum(H(:, k).^2.* P(:, k)) - sum(bag{bagiter,k}*(H(:, k).^2.* P(:, k)))))
              
        % Solve Convex Power Allocation Subproblem (EDM Step)
        cvx_begin quiet
            variable P_new(L, K) nonnegative;
            %maximize sum(W *  log(1 + sum(P_new .* H.^2, 1))/log(2)) - lambda * (sum(P_new(:)) + sum(P_circuit) + sum(P_backhaul));
            maximize  W*sum(1+log(sum(bag{bagiter,k}*(H(:, k).^2.* P_new(:, k)))))/log(2) + W*sum(log(sigma2 + sum(H(:, K).^2.* P_new(:, K)) - sum(bag{bagiter,K}*(H(:, K).^2.* P_new(:, K))))/log(2))  - lambda * (sum(P_new(:)) + sum(P_circuit) + sum(P_backhaul));
            subject to
                %0.1*P_max <= P_new <= P_max;
                %1.1  <= sum(P_new, 2) <= P_max;  % Power constraint at each AP
                P_new(L, K) >= 1;
                sum(P_new, 2) <= P_max;
                %W * log(1 + P_new .* H.^2) >= 1000;
                %W * log(1 + SINR_new(K))/log(2) >= 1000;
        cvx_end
        % Update Power Allocation
        P = P_new;
    end
    
    %% Display Results
    disp(['Bag: ', num2str(bagiter), '']);
    disp('Cluster:');
    disp(bag(bagiter, :));
    disp(['Optimized EE: ', num2str(EE,'%.2e'), ' bits/Joule']);
    disp(['Total Power Consumption: ', num2str(P_total), ' W']);
    disp(['Total Sum-Rate: ', num2str(R_sum / 1e6), ' Mbps']);
    disp('Power Matrix:');
    disp(P);
    disp(['Num Iter: ', num2str(iter),' iterations']);
    disp('');
end