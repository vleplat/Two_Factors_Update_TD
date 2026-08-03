function [X, Z, F1, F2, S, loss_history, primal_residual_history, T, rho] = linreg_krp(Ymn, Phi, rho, mu, szXnm, X0, maxiters, T, verbose)
    % ADMM for
    %   min_Z 0.5*||Ymn - Phi*Z'||_F^2 + 0.5*mu*||Z||_F^2
    %   subject to X = Z and X Khatri-Rao structured.
    %
    % rho is the augmented-Lagrangian penalty beta = 1/gamma used in the
    % manuscript. The update order is X (nonsmooth), Z (smooth), then T.

    if nargin < 9 || isempty(verbose)
        verbose = false;
    end
    if nargin < 8 || isempty(T)
        T = zeros(szXnm(1) * szXnm(2), szXnm(3));
    end

    R = szXnm(3);
    gram = Phi' * Phi + (rho + mu) * eye(R);
    data_term = Ymn' * Phi;

    normY2 = max(norm(Ymn, 'fro')^2, eps);

    X = X0;
    Z = X0;
    loss_history = zeros(maxiters, 1);
    primal_residual_history = zeros(maxiters, 1);

    F1 = zeros(szXnm(1), R);
    F2 = zeros(szXnm(2), R);
    S = ones(R, 1);

    for kiter = 1:maxiters
        Z_old = Z;

        % 1) Update X: exact projection onto the rank-at-most-one set.
        D = Z - T;
        for r = 1:R
            Hr = reshape(D(:, r), szXnm(1), szXnm(2));
            [Ur, Sr, Vr] = svd(Hr, 'econ');
            sigma = max(Sr(1, 1), 0);
            root_sigma = sqrt(sigma);

            % Balanced factor recovery, as stated in the revised paper.
            F1(:, r) = root_sigma * Ur(:, 1);
            F2(:, r) = root_sigma * Vr(:, 1);
            X(:, r) = reshape(F1(:, r) * F2(:, r)', [], 1);
        end

        % 2) Update Z: smooth quadratic subproblem.
        Z = (data_term + rho * (X + T)) / gram;

        % 3) Update the scaled dual variable T.
        primal_residual = X - Z;
        T = T + primal_residual;

        % Relative squared fitting error used by the original scripts.
        loss = norm(Ymn - Phi * Z', 'fro')^2 / normY2;
        dual_residual = rho * (Z - Z_old); %#ok<NASGU>

        loss_history(kiter) = loss;
        primal_residual_history(kiter) = ...
            norm(primal_residual, 'fro') / max(norm(Z, 'fro'), eps);

        if verbose
            fprintf('kiter %d | f = %.6e | d(Z,X) %.6e | rho %.6e\n', ...
                kiter, loss, primal_residual_history(kiter), rho);
        end

        if (kiter > 5) && (loss <= 1e-6) && ...
                (primal_residual_history(kiter) <= 1e-5)
            break;
        end
    end

    loss_history = loss_history(1:kiter);
    primal_residual_history = primal_residual_history(1:kiter);
end
