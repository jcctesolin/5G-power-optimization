function [channelofCFsetups, largeScaleFading1, channelofCFsetups_tilde]=channelcellfree_N(APcellfree,UElocations,nbrOfAPs,N,K,pack)
% This function generates cell-free channel for nbrOfAPs as number of APs,
% each AP equipped with N antennas, and K is number of UEs

channelofCFsetups1 = zeros(nbrOfAPs*N,K); 
channelofCFsetups_tilde= zeros(nbrOfAPs*N,K);

SNR = @(hor_dist) db2pow(10+96-30.5-36.7*log10(sqrt(hor_dist.^2+10^2)));
largeScaleFading = zeros(nbrOfAPs, K); 
for k = 1:K
    for n = 1:N
        % Compute distance for each AP to the UE, considering the APs having multiple antennas
        % The antennas at the same AP have the same distance to the UE
        distanceCellfree = abs(APcellfree(:) - UElocations(pack,k));
        
        % Apply large scale fading, which is the same for all antennas at each AP
        % largeScaleFading = sqrt(SNR(distanceCellfree));
        largeScaleFading(:, k) = SNR(distanceCellfree);
        % Generate random small scale fading for each antenna at the AP
        smallScaleFading = exp(1i*2*pi*rand(nbrOfAPs,1));
        smallScaleFading_tilde = exp(1i*2*pi*rand(nbrOfAPs,1));
        
        % Index to represent the antenna in the channel matrix
        idx = (n-1)*nbrOfAPs + (1:nbrOfAPs);
        
        channelofCFsetups1(idx,k) = sqrt(0.95)*sqrt(largeScaleFading(:, k)) .*smallScaleFading;
        channelofCFsetups_tilde(idx,k) = sqrt(0.05)*sqrt(largeScaleFading(:, k)) .*smallScaleFading_tilde;
    end
end

channelofCFsetups=transpose(channelofCFsetups1+channelofCFsetups_tilde); %Channel cell-free including both estimated and estimation error parts
channelofCFsetups_tilde=transpose(channelofCFsetups_tilde); % Channel estimation error 
largeScaleFading1=transpose(largeScaleFading);