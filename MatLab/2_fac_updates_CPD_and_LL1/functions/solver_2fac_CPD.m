function [Y_hat, mainloss_history] = solver_2fac_CPD(Y, R, Y_hat, rho, mu, maxoutiters, maxiters, min_rho_stable)
    % Two-factor CPD solver.
    %
    % The inner ADMM uses the order X (nonsmooth), Z (smooth), then T.
    % NOTE: this routine keeps the submitted numerical schedule of overlapping
    % consecutive pairs (1,2), (2,3), ... . This differs from the random,
    % non-overlapping pairing shown in Algorithm 2 and is left unchanged here
    % to preserve the numerical experiment design.

    modes = ndims(Y);
    szY = size(Y);
    normY = max(frob(Y), eps);
    normY2 = normY^2;
    mainloss_history = zeros(maxoutiters * max(modes - 1, 1), 1);
    history_idx = 0;

    rho_stable = rho;
    counter = 0;
    flag = 0;

    for kiter = 1:maxoutiters
        for n = 1:modes-1
            m = n + 1;

            modes_1 = sort([n, m]);
            modes_2 = sort(setdiff(1:modes, modes_1));

            sz2 = prod(szY(modes_2));
            sz1 = prod(szY(modes_1));

            Y_nm = permute(Y, [modes_2, modes_1]);
            Y_nm = reshape(Y_nm, [sz2, sz1]);

            szXnm = [szY(modes_1(1)), szY(modes_1(2)), R];

            U_in = Y_hat.factors{modes_1(1)};
            V_in = Y_hat.factors{modes_1(2)};
            X0 = kr(V_in, U_in * diag(Y_hat.weights));

            Factorx2 = Y_hat.factors(modes_2(end:-1:1));
            Phi = kr(Factorx2);

            fcurr = norm(Y_nm - Phi * X0', 'fro')^2 / normY2;

            best_fit = inf;
            accepted = false;
            max_retries = 20;
            for retry = 1:max_retries
                [~, ~, Unew, Vnew, Snew, loss_history, dZ, T, ~] = ...
                    linreg_krp(Y_nm, Phi, rho, mu, szXnm, X0, ...
                    maxiters, [], false); %#ok<ASGLU>

                X_candidate = kr(Vnew, Unew * diag(Snew));
                fit_candidate = norm(Y_nm - Phi * X_candidate', 'fro')^2 / normY2;

                if fit_candidate < best_fit
                    best_fit = fit_candidate;
                    best_Unew = Unew;
                    best_Vnew = Vnew;
                    best_Snew = Snew;
                    best_loss_history = loss_history;
                    best_dZ = dZ;
                end

                if fit_candidate <= max(10 * fcurr, fcurr + 1e-14)
                    accepted = true;
                    break;
                end
                rho = 1.1 * rho;
            end

            if ~accepted
                warning('solver_2fac_CPD:RetryLimit', ...
                    ['The pair update did not pass the acceptance test after %d retries. ' ...
                     'The best trial is retained.'], max_retries);
                Unew = best_Unew;
                Vnew = best_Vnew;
                Snew = best_Snew;
                loss_history = best_loss_history;
                dZ = best_dZ;
            end

            Y_hat.factors{modes_1(1)} = Unew;
            Y_hat.factors{modes_1(2)} = Vnew;
            Y_hat.weights = Snew;

            % Record the actual relative CP reconstruction error associated
            % with the structured factors X, rather than the auxiliary Z-fit.
            factors_eval = Y_hat.factors;
            factors_eval{1} = factors_eval{1} * diag(Y_hat.weights);
            current_relerr = frob(Y - cpdgen(factors_eval)) / normY;

            history_idx = history_idx + 1;
            mainloss_history(history_idx) = current_relerr;

            fprintf('kiter %d | pair (%d,%d) | relerr = %.6e | d(Z,X) %.6e | rho %.6e\n', ...
                kiter, n, m, current_relerr, dZ(end), rho);

            if (history_idx > 5) && ...
                    (current_relerr <= min(mainloss_history(history_idx-5:history_idx-1)))
                counter = counter + 1;
                flag = 1;
            else
                counter = 0;
                flag = 0;
            end

            if (counter == 2) && (dZ(end) < 1e-6)
                rho = 0.95 * rho;
                counter = 0;
            end

            if (history_idx > 1) && ...
                    (current_relerr > mainloss_history(history_idx-1))
                counter = counter - 1;
            end

            if counter < 0
                rho = min(rho + 0.1, rho_stable);
                counter = 0;
            end

            if flag == 1
                rho_stable = max(min_rho_stable, min(rho, rho_stable));
            end
        end

        if (kiter > 5) && (mainloss_history(history_idx) <= 1e-6) && ...
                (dZ(end) <= 1e-4)
            break;
        end
    end

    mainloss_history = mainloss_history(1:history_idx);
end
