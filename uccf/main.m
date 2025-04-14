clc; clear;

% Parâmetros do sistema
K = 2;      % número de usuários
L = 2;      % número de APs
sigma2 = 1e-9; % ruído
P_max = 0.2 * ones(L,1);  % potência máxima por AP (W)
Pc = 0.1 * L;             % consumo fixo dos circuitos (W)

% Ganhos de canal aleatórios h(k,l)
h = abs(randn(K,L) + 1j*randn(K,L)).^2;

% Inicialização
p0 = 0.1 * ones(L,1);   % palpite inicial
lambda = 0;             % estimativa inicial da eficiência
tol = 1e-4;
max_iter = 50;
iter = 0;
eta_diff = Inf;

options = optimoptions('fmincon','Display','none','Algorithm','sqp');

while eta_diff > tol && iter < max_iter
    obj = @(p) -sum_rate(p,h,sigma2,K,L) + lambda * total_power(p,Pc);
    
    [p_opt, fval] = fmincon(obj, p0, [], [], [], [], zeros(L,1), P_max, [], options);
    
    R = sum_rate(p_opt,h,sigma2,K,L);
    P = total_power(p_opt,Pc);
    eta_new = R / P;
    
    eta_diff = abs(eta_new - lambda);
    lambda = eta_new;
    p0 = p_opt;
    iter = iter + 1;
    
    fprintf('Iteração %d: EE = %.4f bit/J\n', iter, eta_new);
end

fprintf('\nSolução final:\n');
disp('Potência ótima dos APs:');
disp(p_opt);
fprintf('Eficiência Energética final: %.4f bit/J\n', lambda);
