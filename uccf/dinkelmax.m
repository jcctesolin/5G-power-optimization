function [x_opt, f_opt, iter] = dinkelmax(f, g, x0, lb, ub, tol, maxIter)
    % 
    % Dinkelbach algorithm for fractional optimization
    %
    % Input:
    %   f       - numerator f(x)
    %   g       - denominator g(x)
    %   x0      - initial guess
    %   lb, ub  - lower and upper variable boundaries
    %   tol     - convergence tolerance
    %   maxIter - max iterations
    %
    % Output:
    %   x_opt   - optimal solution 
    %   f_opt   - optimal value of  f(x)/g(x)
    %   iter    - num of iterations performed

    % initial λ (lambda) as f(x0) / g(x0)
    lambda = f(x0) / g(x0);
    iter = 0;
    options = optimset('Display', 'off'); % Option of silent optimization

    while iter < maxIter
        iter = iter + 1;

        % solve optimization problem to find x*
        objFunc = @(x) -(f(x) - lambda * g(x)); 
        x_opt = fmincon(objFunc, x0, [], [], [], [], lb, ub, [], options);

        % update lambda
        f_x = f(x_opt);
        g_x = g(x_opt);
        new_lambda = f_x / g_x;

        % stop criteria
        if abs(new_lambda - lambda) < tol
            break;
        end

        lambda = new_lambda; % update lambda for the next iteration
    end

    f_opt = f(x_opt) / g(x_opt); % fractional function optimal value
end
