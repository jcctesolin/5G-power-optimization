clc; clear; close all;

% Parâmetros do sistema
L = 4;  % Número total de APs
K = 2;  % Número total de usuários
P_max = 10;  % Potência máxima por AP
Pc = 2; % Consumo fixo de energia
B = 1e3; % Largura de banda (Hz)
R_min = 100 / B; % Throughput mínimo exigido (bps transformado para bps/Hz)
tol = 1e-6; % Critério de parada
lambda = 0; % Inicialização do parâmetro de Dinkelbach
max_iter = 1000; % Limite máximo de iterações
max_srv = 3; % Max number of serving APs

% Matriz de canais Rayleigh para todos os APs e usuários
H_full = abs(randn(L, K));  
% Inicialização da potência dos APs
P = ones(L,K) * (P_max / 2);

%% Creating clusters 
D = zeros(L);

[sorted_vals, sorted_indices] = sort(H_full);

for j=1:max_srv
    for i=1:K
        if j = 1 
            then mst_AP_index = sorted_indices(j:i);
        D(sorted_indices(j:i),K) = 1;
        cluster=[sorted_indices(j:i)]; 
    end
end  

% Criar Clusters de APs com diferentes tamanhos (exemplo)
clusters = { [1], [1,2], [2,3,4], [1,2,3,4,5,6] }; % Diferentes clusters de APs




%%
% Função para calcular SINR por usuário considerando apenas APs do cluster dele
SINR = @(P) arrayfun(@(k) ...
    sum(H_full(k, clusters{k}) .* P(clusters{k})) / ...
    (0.1 + sum(H_full(k, :) .* P(:)) - sum(H_full(k, clusters{k}) .* P(clusters{k}))), ...
    1:K);

% Função de taxa de transmissão (R)
R = @(P) sum(log2(1 + SINR(P)));

% Função de consumo de energia (E)
E = @(P) sum(P) %+ Pc;

% Algoritmo de Dinkelbach
iter = 0;
while true
    iter = iter + 1;

    % Resolver subproblema: maximizar R(P) - lambda * E(P) com restrições
    cvx_begin quiet
        variable P(L,1) nonnegative;
        maximize( sum(log2(1 + arrayfun(@(k) ...
            sum(H_full(k, clusters{k}) .* P(clusters{k})) / ...
            (0.1 + sum(H_full(k, :) .* P(:)) - sum(H_full(k, clusters{k}) .* P(clusters{k}))), ...
            1:K))) - lambda * (sum(P) + Pc) )
        subject to
            0 <= P <= P_max; % Restrição de potência por AP
            log2(1 + SINR(P)) >= R_min; % Restrição de throughput mínimo
    cvx_end

    % Atualizar valores
    R_val = R(P);
    E_val = E(P);
    lambda_new = R_val / E_val;

    % Critério de parada
    if abs(R_val - lambda * E_val) < tol || iter > max_iter
        break;
    end
    lambda = lambda_new;
end

% Determinar o cluster final associado a cada usuário
final_clusters = cell(K,1);
for k = 1:K
    final_clusters{k} = clusters{k}(P(clusters{k}) > 0); % Apenas APs com potência alocada
end

% Exibir resultados
disp('Potência ótima dos APs:');
disp(P);
disp(['Eficiência energética ótima: ', num2str(R(P) / E(P))]);
disp(['Número de iterações: ', num2str(iter)]);

disp('Clusters finais dos usuários:');
for k = 1:K
    fprintf('Usuário %d: APs ativos -> [%s]\n', k, num2str(final_clusters{k}));
end

