function R = sum_rate_mrt(p, h, noise_power)
[K, L, M] = size(h);
R = 0;
for k = 1:K
    signal = 0;
    interf = 0;
    for l = 1:L
        hk = squeeze(h(k, l, :));
        signal = signal + sqrt(p(k,l)) * norm(hk);
    end
    for j = 1:K
        if j ~= k
            for l = 1:L
                hj = squeeze(h(j, l, :));
                interf = interf + p(j,l) * norm(hj)^2;
            end
        end
    end
    SINR_k = (signal^2) / (interf + noise_power);
    R = R + log(1 + SINR_k);
end