function [Y_hat, mainloss_history] = solver_2fac_CPD(Y, R, Y_hat, rho, mu, maxoutiters, maxiters, min_rho_stable)
    % maxoutiters = 100
    % rho = 1e1
    % maxiters = 10
    modes = ndims(Y);
    szY = size(Y); 

    normY = norm(reshape(Y, [], 1))^2;
    mainloss_history = [];

    rho_stable = rho;
    counter = 0;
    flag = 0;

    for kiter = 1:maxoutiters
        for n = 1:modes-1
            m = n + 1; % can change the mode m 

            modes_1 = sort([n, m]);
            modes_2 = sort(setdiff(1:modes, modes_1));

            sz2 = prod(szY(modes_2));
            sz1 = prod(szY(modes_1));

            % unfolding mode-(n,m)    
            Y_nm = permute(Y, [modes_2, modes_1]);
            Y_nm = reshape(Y_nm, [sz2, sz1]);

            szXnm = [szY(modes_1(1)), szY(modes_1(2)), R];

            % initialize X0
            U_in = Y_hat.factors{modes_1(1)};
            V_in = Y_hat.factors{modes_1(2)};
            X0 = kr(V_in, U_in * diag(Y_hat.weights)); % Replace with your kr function

            % Compute Phi
            Factorx2 = Y_hat.factors(modes_2(end:-1:1));
            Phi = kr(Factorx2); % Replace with your kr function

            % solve sub-problem
            fcurr = norm(Y_nm - Phi * X0', 'fro')^2 / normY;

            while true
                [X, Z, Unew, Vnew, Snew, loss_history, dZ, T, rhonew] = linreg_krp(Y_nm, Phi, rho, mu, szXnm, X0, maxiters, [], false);

                if loss_history(end) > fcurr * 10
                    % fail case
                    rho = rho * 1.1;
                else
                    break;
                end
            end

            % Update factors in Y_hat
            Y_hat.factors{modes_1(1)} = Unew;
            Y_hat.factors{modes_1(2)} = Vnew;
            Y_hat.weights = Snew;                

            fprintf('kiter %2.2f | f = %2.4e | d(Z,X) %2.4e | rho % .4e\n', kiter, loss_history(end), dZ(end), rho_stable);

            mainloss_history = [mainloss_history; loss_history(end)];

            % If the loss function value is smaller than the previous 10 values, increment the counter by one and set the flag to True
            if (length(mainloss_history) > 5) && (kiter > 1) && (loss_history(end) <= min(mainloss_history(end-5:end)))
                counter = counter + 1;
                flag = 1;
            else
                counter = 0;
                flag = 0;
            end

            % If the counter reaches 10 and the error is smaller than 1e-6, reduce rho by a factor of 0.9 or 0.95, and reset the counter to zero
            if (counter == 2) && (dZ(end) < 1e-6)
                rho = rho * 0.95;
                counter = 0;
            end

            % If the loss function value is larger than the previous value, reduce the counter by one
            if (length(mainloss_history) > 2) && (kiter > 1) && (loss_history(end) > mainloss_history(end-1))
                counter = counter - 1;
            end

            % If the counter is negative, increase rho by a small amount (e.g., 0.1) until it reaches rho_stable, and reset the counter to zero
            if counter < 0
                rho = min(rho + 0.1, rho_stable);
                counter = 0;
            end

            % Update rho_stable to be the minimum of rho and rho_stable only when the flag is True
            if flag == 1
                rho_stable = max(min_rho_stable, min(rho, rho_stable));
            end
        end

        if (kiter > 5) && (mainloss_history(end) <= 1e-6) && (dZ(end) <= 1e-4)
            break;
        end
    end
end
