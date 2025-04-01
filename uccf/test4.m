clc; clear; close all;

% Parâmetros do sistema
M = 6;  % Número de APs
K = 4;  % Número de usuários
P_max = 10;  % Potência máxima por AP
Pc = 2; % Consumo fixo de energia
tol = 1e-6; % Critério de parada
lambda = 0; % Inicialização do parâmetro de Dinkelbach

% Ganhos de canal entre APs e usuários (modelo Rayleigh)
H = abs(randn(K, M));  

% Ruído branco gaussiano
sigma2 = 0.1;  

% Inicialização da potência dos APs
P = ones(M,1) * (P_max / 2);

% Função SINR considerando interferência
SINR = @(P) diag(H * P) ./ (sigma2 + sum(H * P, 2) - diag(H * P));

% Função de taxa de transmissão (R)
R = @(P) sum(log2(1 + SINR(P)));

% Função de consumo de energia (E)
E = @(P) sum(P) + Pc;

% Algoritmo de Dinkelbach para Cell-Free
iter = 0;
max_iter = 1000; 

while true
    iter = iter + 1;

    % Resolver subproblema: max R(P) - lambda * E(P)
    cvx_begin quiet
        variable P(M,1) nonnegative;
        maximize( sum(log2(1 + diag(H * P) ./ (sigma2 + sum(H * P, 2) - diag(H * P)))) - lambda * (sum(P) + Pc) )
        subject to
            0 <= P <= P_max; 
    cvx_end

    % Atualizar valores
    R_val = R(P);
    E_val = E(P);

    % Atualizar lambda
    lambda_new = R_val / E_val;

    % Critério de parada
    if abs(R_val - lambda * E_val) < tol || iter > max_iter
        break;
    end

    lambda = lambda_new;
end

% Exibir resultados
disp('Potência ótima dos APs:');
disp(P);
disp(['Eficiência energética ótima: ', num2str(R(P) / E(P))]);
disp(['Número de iterações: ', num2str(iter)]);











