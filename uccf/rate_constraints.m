function [c, ceq] = rate_constraints(p, h, noise_power, Rmin)
[K, L, M] = size(h);
c = zeros(K,1);
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
    rate_k = log(1 + SINR_k);
    c(k) = Rmin - rate_k;
end
ceq = [];