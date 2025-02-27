function A_APs_selected = APselection(largeScaleFading,N)% This function performs AP selection using LSF criterion
    
    % Get large scale fading and channel coefficients from LSFforcellfree
    
    Magn = largeScaleFading;  % Replaced LSF with largeScaleFading
    

%     % Assuming you already have largeScaleFading matrix defined
[K, L] = size(largeScaleFading);  % Get dimensions of the matrix
% [K, M] = size(P1mmse);
a = zeros(K, L*N);  % Pre-allocation for efficiency

alpha_lsf = sum(sum(largeScaleFading)) / (L*N*K);


for ii = 1:K
    for jj = 1:L
        if Magn(ii, jj) >= alpha_lsf
            % If the condition is met, set the corresponding N elements to 1
            a(ii, (jj-1)*N+1:jj*N) = 1;
        end
    end
    A(:,:,ii) = diag(a(ii,:));
end
    
 A_APs_selected = A;