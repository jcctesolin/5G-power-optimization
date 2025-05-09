% Arquivo: gerar_channel_data_mimo.m
% Este script gera dados de canal MIMO simulados e salva em um arquivo .mat
% Os dados representam a matriz de canais entre usuarios e APs (K x L).

% Parâmetros do sistema
nUE = 2;      % Número de usuários
nAP = 4;      % Número de Access Points
nAntpUE = 1; % Número de antenas por UE
nAntpAP = 4; % Número de antenas por AP

% Inicializar matriz de canais H
% Dimensão: (antennas_per_AP * num_APs) x (antennas_per_UE * num_UEs)
H = (randn(nAntpUE*nUE, nAntpAP*nAP) + 1i*randn(nAntpUE*nUE, nAntpAP*nAP)) / sqrt(2);

% Inserir desvanecimento de larga escala (pathloss)
% Supomos distância aleatória entre 10m e 100m
distances = 10 + 90*rand(nUE, nAP); % Distância UE-AP (em metros)
pathloss_exponent = 3.5; % Expoente típico de pathloss
pathloss = (distances.^(-pathloss_exponent/2)); % Fator de atenuação

% Aplicar pathloss na matriz H
for ue = 1:nUE
    for ap = 1:nAP
        idx_ue = (ue-1)*nAntpUE + (1:nAntpUE);
        idx_ap = (ap-1)*nAntpAP + (1:nAntpAP);
        H(idx_ue, idx_ap) = H(idx_ue, idx_ap) * pathloss(ue, ap);
    end
end

% Salvar a matriz H em arquivo .mat
save('channel_data_mimo2.mat', 'H', 'nUE', 'nAP', 'nAntpUE', 'nAntpAP', 'distances');

disp('Arquivo channel_data_mimo2.mat gerado com sucesso.');
