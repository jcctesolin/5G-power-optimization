function R = sum_rate(p, h, sigma2, K, L)
    R = 0;
    for k = 1:K
        signal = 0;
        interference = 0;
        for l = 1:L
            signal = signal + sqrt(p(l)) * h(k,l);
        end
        for j = 1:L
            interference = interference + (j ~= l) * p(j) * h(k,j);
        end
        SINR_k = abs(signal)^2 / (interference + sigma2);
        R = R + log2(1 + SINR_k);
    end
end
