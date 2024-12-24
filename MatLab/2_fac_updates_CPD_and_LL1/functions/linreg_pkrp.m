function [U, loss_history, primal_residual_history, T, rho] = linreg_pkrp(Ymn, U, L, Phi, rho, mu, szXnm, X0, maxiters, T, verbose)
    % ADMM
    %   solve min||Y-Phi*Z^T||^2_f  
    %   s.t. Z= X 

    % szY(1) = [I1 I2 R]
    R = length(L);
    Q = inv(Phi' * Phi + (rho + mu) * eye(R));
    B = Ymn' * Phi;

    if nargin < 10 || isempty(T)
        T = zeros(szXnm(1) * szXnm(2), R);
    end

    normY = norm(Ymn, 'fro')^2;

    X = X0; 
    Z = X;
    f0 = norm(Ymn - Phi * X', 'fro')^2 / normY;
    loss_history = []; 
    primal_residual_history = [];

    % Define the parameters of the rule
    mu = 10; % Residual ratio threshold
    tau = 2; % Rho scaling factor

    for kiter = 1:maxiters
        Z_old = Z;

        % Update Z
        Z = (B + rho * (X + T)) * Q;

        % Update X
       
        D = Z - T;
        F1 = zeros(szXnm(1), R);
        % D = reshape(Z - T, [szXnm(2), szXnm(1), szXnm(3)]);
        F2 = zeros(szXnm(2), R); 
        S = zeros(R, 1);
        for r = 1:R
            % [u, s, v] = svds(D(:, :, r)', 1);
            Hr = reshape(D(:,r), szXnm(1), szXnm(2));
            [u, s, v] = svds(Hr, L(r));
            % F1(:, r) = u;
            % F2(:, r) = v';
            % S(r) = s;
            % prod = u*s*v';
            % X(:,r) = prod(:);
            U{r}{1} = u*s;
            U{r}{2} = v;
            % U{r}{4} = s;
        end
        
        % X = kr(F2, F1 * diag(S)); % Replace with your kr function
        X = pw_vecL(U,R,L);

        % Compute the primal residual
        r = X - Z;

        % Update T
        T = T + r;

        % Compute the objective function value
        loss = norm(Ymn - Phi * Z', 'fro')^2 / normY; % f(Z) + i(X) 

        % dual residual
        s = rho * (Z - Z_old);
        % Compute the primal and dual residual norms
        norm_r = norm(r, 'fro');
        norm_s = norm(s, 'fro');  

        loss_history = [loss_history; loss];
        primal_residual_history = [primal_residual_history; norm_r / norm(Z, 'fro')];

        if verbose
            fprintf('kiter %d | f = %f | d(Z,X) %f | rho %f\n', kiter, loss, primal_residual_history(end), rho);
        end 

        if (kiter > 5) && (loss <= 1e-6) && (primal_residual_history(end) <= 1e-5)
            break; 
        end 
    end
end