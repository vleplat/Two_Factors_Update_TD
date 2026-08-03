function [Y_hat, U, mainloss_history, U0] = solver_2fac_ll1_zfirst_reference(Y, L, U, rho, mu, maxoutiters, maxiters, min_rho_stable, init_type, varargin)
    % Two-factor update solver for an (L_r,L_r,1)-BTD.
    %
    % The A,B inner update uses mu_AB = 0, consistently with the revised
    % manuscript. The optional ridge parameter mu is used only in the C-update.

    if nargin < 3
        U = [];
    end
    if nargin < 4 || isempty(rho)
        rho = 10; % augmented-Lagrangian penalty beta = 1/gamma
    end
    if nargin < 5 || isempty(mu)
        mu = 0;
    end
    if nargin < 6 || isempty(maxoutiters)
        maxoutiters = 200;
    end
    if nargin < 7 || isempty(maxiters)
        maxiters = 100;
    end
    if nargin < 8 || isempty(min_rho_stable)
        min_rho_stable = 0.1;
    end
    if nargin < 9 || isempty(init_type)
        init_type = 'RAND';
    end

    szY = size(Y);
    normY = max(frob(Y), eps);
    normY2 = normY^2;
    mainloss_history = zeros(maxoutiters, 1);

    rho_stable = rho;
    counter = 0;
    flag = 0;
    R = length(L);

    % Check/format the tensor for Tensorlab routines.
    type = getstructure(Y);
    isstructured = ~any(strcmp(type, {'full', 'sparse', 'incomplete'}));
    if ~isstructured
        Y = fmt(Y);
    end
    size_tens = getsize(Y);

    % Initialize the factor matrices unless provided by the user.
    fprintf('Step 2: Initialization \n');
    if isempty(U)
        if strcmpi(init_type, 'RAND')
            fprintf('is ll1_rnd (default) \n');
            [U, ~] = ll1_rnd(size_tens, L, varargin{:});
        elseif strcmpi(init_type, 'GEVD')
            fprintf('is ll1_gevd \n');
            [U, ~] = ll1_gevd(Y, L, varargin{:});
        else
            error('solver_2fac_ll1_zfirst_reference:UnknownInitialization', ...
                'Unknown initialization type "%s".', init_type);
        end
    else
        fprintf('is manual... \n');
    end
    U0 = U;

    % Mode-3 unfolding.
    idx = 1:3;
    mode = 3;
    Y_3 = tens2mat(Y, mode, idx(idx ~= mode));
    C = zeros(szY(3), R);
    for r = 1:R
        C(:, r) = U{r}{3};
    end

    for kiter = 1:maxoutiters
        % ---------------------------------------------------------------
        % Update A=[A_1,...,A_R] and B=[B_1,...,B_R].
        % The revised derivation uses no ridge term on the merged variable.
        % ---------------------------------------------------------------
        X0 = pw_vecL(U, R, L);
        fcurr = norm(Y_3 - C * X0', 'fro')^2 / normY2;

        U_before = U;
        best_fit = inf;
        best_U = U;
        best_loss_history = [];
        best_dZ = [];
        best_T = [];
        accepted = false;
        max_retries = 20;

        for retry = 1:max_retries
            U_trial = U_before;
            [U_trial, loss_history, dZ, T, ~] = ...
                linreg_pkrp_zfirst_reference(Y_3, U_trial, L, C, rho, 0, ...
                getsize(Y), X0, maxiters, [], false);

            X_trial = pw_vecL(U_trial, R, L);
            fit_trial = norm(Y_3 - C * X_trial', 'fro')^2 / normY2;

            if fit_trial < best_fit
                best_fit = fit_trial;
                best_U = U_trial;
                best_loss_history = loss_history;
                best_dZ = dZ;
                best_T = T;
            end

            if fit_trial <= max(10 * fcurr, fcurr + 1e-14)
                U = U_trial;
                accepted = true;
                break;
            end

            rho = 1.1 * rho;
        end

        if ~accepted
            warning('solver_2fac_ll1_zfirst_reference:RetryLimit', ...
                ['The A,B inner update did not pass the acceptance test after %d retries. ' ...
                 'The best trial is retained.'], max_retries);
            U = best_U;
            loss_history = best_loss_history;
            dZ = best_dZ;
            T = best_T; %#ok<NASGU>
        end

        % Update C in ALS/ridge-regression fashion.
        Phi = pw_vecL(U, R, L);
        gramC = Phi' * Phi + mu * eye(R);
        C = (gramC \ (Y_3 * Phi)')';
        for r = 1:R
            U{r}{3} = C(:, r);
        end

        current_relerr = frob(btdres(Y, U)) / normY;
        mainloss_history(kiter) = current_relerr;

        fprintf('kiter %d | relerr = %.6e | d(Z,X) %.6e | rho %.6e\n', ...
            kiter, current_relerr, dZ(end), rho);

        % Keep the original empirical outer penalty strategy, but compare
        % quantities with consistent units (relative error vs relative error).
        if (kiter > 5) && ...
                (current_relerr <= min(mainloss_history(kiter-5:kiter-1)))
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

        if (kiter > 1) && (current_relerr > mainloss_history(kiter-1))
            counter = counter - 1;
        end

        if counter < 0
            rho = min(rho + 0.1, rho_stable);
            counter = 0;
        end

        if flag == 1
            rho_stable = max(min_rho_stable, min(rho, rho_stable));
        end

        if (kiter > 5) && (current_relerr <= 1e-6) && (dZ(end) <= 1e-4)
            break;
        end
    end

    mainloss_history = mainloss_history(1:kiter);
    Y_hat = ful(U);
end
