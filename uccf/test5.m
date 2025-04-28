L = 4; % Tamanho do vetor desejado

% Gera todas as combinações possíveis de 0 e 1, excluindo o vetor de zeros
matriz_combinacoes = dec2bin(1:(2^L - 1)) - '0';

% Converte cada linha em um vetor separado dentro de uma célula
vetores = mat2cell(matriz_combinacoes, ones(size(matriz_combinacoes,1),1), L);

% Exibe os vetores
for i = 1:length(vetores)
    disp(vetores{i})
end

