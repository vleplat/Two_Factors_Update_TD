function [U, loss_history, primal_residual_history, T, rho] = linreg_pkrp_zfirst_reference(Ymn, U, L, Phi, rho, mu, szXnm, X0, maxiters, T, verbose)
    % Reference implementation used only to compare update orders.
    % It is identical to the audited X-first routine except that it uses the
    % former order Z (smooth), X (nonsmooth), then T.

    if nargin < 11 || isempty(verbose)
        verbose = false;
    end
    if nargin < 10 || isempty(T)
        T = zeros(szXnm(1) * szXnm(2), length(L));
    end

    R = length(L);
    gram = Phi' * Phi + (rho + mu) * eye(R);
    data_term = Ymn' * Phi;
    normY2 = max(norm(Ymn, 'fro')^2, eps);

    X = X0;
    Z = X0;
    loss_history = zeros(maxiters, 1);
    primal_residual_history = zeros(maxiters, 1);

    for kiter = 1:maxiters
        Z_old = Z;

        % Former order: Z first.
        Z = (data_term + rho * (X + T)) / gram;

        % Then X.
        D = Z - T;
        for r = 1:R
            Hr = reshape(D(:, r), szXnm(1), szXnm(2));
            if L(r) > min(size(Hr))
                error('linreg_pkrp_zfirst_reference:InvalidRank', ...
                    'L(%d)=%d exceeds min(size(H_r))=%d.', ...
                    r, L(r), min(size(Hr)));
            end
            [Ur, Sr, Vr] = svd(Hr, 'econ');
            Ur = Ur(:, 1:L(r));
            Vr = Vr(:, 1:L(r));
            singular_values = max(diag(Sr(1:L(r), 1:L(r))), 0);
            root_Sr = diag(sqrt(singular_values));
            U{r}{1} = Ur * root_Sr;
            U{r}{2} = Vr * root_Sr;
        end
        X = pw_vecL(U, R, L);

        primal_residual = X - Z;
        T = T + primal_residual;

        loss = norm(Ymn - Phi * Z', 'fro')^2 / normY2;
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
