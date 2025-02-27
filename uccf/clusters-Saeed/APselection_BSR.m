function A_APs_selected = APselection_BSR(channelCellfree, channelCellfree_tilde, Es, sigma)% This function performs AP selection using BSR criterion
    % Compute MMSE beamforming coefficients
    [~, Pmmse] = beamformingMMSE_I(channelCellfree, channelCellfree_tilde, Es, sigma);

    % Calculate the direct channel
    channelCellfree_hat = channelCellfree - channelCellfree_tilde;

    % Extract dimensions from the Pmmse matrix
    [M, K] = size(Pmmse);

    % Pre-allocate SR matrix
    SR = zeros(K, M);

    S1 = 0;
    for k = 1:K
        S2 = 0;
        for m = 1:M
            gkmh = abs(channelCellfree_hat(k, m));
            gkmt = abs(channelCellfree_tilde(k, m));
            Pk = Pmmse(:, k);
            power = Pk' * Pk; % Power due to precoder
            snr = (gkmh^2 * power) / (gkmt^2 * power + sigma^2);
            SR(k, m) = log2(1 + snr);
            S2 = S2 + SR(k, m);
        end
        S1 = S1 + S2;
    end

    alpha_src = 1/(K * M) * S1;

    % Step 1: ASR Calculation and Initial AP Selection
    a = zeros(K, M);
    for k = 1:K
        a(k, SR(k, :) >= alpha_src) = 1;
    end
    
    % Step 2: Evaluation of AP Coverage
    N_AP_Av = mean(sum(a, 2));
    
    % Step 3: Identification of Under-supported UEs
    S_Low = find(sum(a, 2) < N_AP_Av);
    
    % Step 4: Augmenting AP Coverage
    for k = S_Low'
        while sum(a(k, :)) < N_AP_Av
            [~, m_max] = max(SR(k, ~a(k, :)));
            m_max_idx = find(~a(k, :));
            a(k, m_max_idx(m_max)) = 1;
        end
    end
    
    % Storing the AP selection matrix
    A_APs_selected = zeros(M, M, K);
    for ii = 1:K
        A_APs_selected(:,:,ii) = diag(a(ii,:));
    end
end