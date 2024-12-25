function [Y_hat, U, mainloss_history, U0] = solver_2fac_ll1(Y, L, U, rho, mu, maxoutiters, maxiters, min_rho_stable, init_type, varargin)

    % STEP0: Check input parameters and misc.
    min_rho_stable = 0.1;

    if ~exist('mu', 'var') || isempty(mu)
        mu = 0;
    end

    if ~exist('maxoutiters', 'var') || isempty(maxoutiters)
        maxoutiters = 200;
    end

    if ~exist('maxiters', 'var') || isempty(maxiters)
        maxiters = 100;
    end

    if ~exist('rho', 'var') || isempty(rho)
        rho = 2; %small:2, mid:10
    end

    if ~exist('init_type', 'var') || isempty(init_type)
        init_type = 'RAND';
    end

    % modes = ndims(Y);
    szY = size(Y); 

    normY = norm(reshape(Y, [], 1))^2;
    mainloss_history = [];

    rho_stable = rho;
    counter = 0;
    % flag = 0;
    
    R = length(L);

    % STEP1: Check the tensor Y.
    type = getstructure(Y);
    isstructured = ~any(strcmp(type, {'full', 'sparse', 'incomplete'}));
    
    if ~isstructured
        Y = fmt(Y);
        type = getstructure(Y);
    end
    size_tens = getsize(Y);

    % STEP 2: initialize the factor matrices unless they were all provided by
    % the user.
    fprintf('Step 2: Initialization \n');
    if ~exist('U', 'var') || isempty(U)
        if strcmp(init_type,'RAND')
            fprintf('is ll1_rnd (default) \n');
            [U,~] = ll1_rnd(size_tens, L, varargin);
        elseif strcmp(init_type,'GEVD')
            % GEVD not working in this current version - TBD
            fprintf('is ll1_gevd \n');
            [U, ~] = ll1_gevd(Y, L, varargin);
        end
    else
        fprintf('is manual... \n');
    end
    U0 = U;
    
    % STEP 3: 3-mode Unfolding
    idx=1:3;
    mode = 3;
    Y_3 = tens2mat(Y,mode,idx(idx~=mode));
    C = zeros(szY(3),R);
    for r=1:R
        C(:,r) = U{r}{3};
    end

    % Main Loop
    for kiter = 1:maxoutiters

        % Update of A=[A_1,...,A_R] and B=[B_1,...,B_R]
        X0 = pw_vecL(U,R,L);
        Phi = C;
        fcurr = norm(Y_3 - Phi * X0', 'fro')^2 / normY;
        while true
            [U, loss_history, dZ, T, rhonew] = linreg_pkrp(Y_3, U, L, C, rho, mu, getsize(Y), X0, maxiters, [], false);

            if loss_history(end) > fcurr * 10
                % fail case
                rho = rho * 1.1;
            else
                break;
            end
        end

        
        % Update of C; ALS fashion
        %%% Compute Phi
        Phi = pw_vecL(U,R,L);
        % C = (Y_3*Phi)/(Phi'*Phi+mu*eye(R));
        Ct = linsolve(Phi'*Phi+mu*eye(R),(Y_3*Phi)');
        C = Ct';
        for r=1:R
            U{r}{3} = C(:,r);
        end

        fprintf('kiter %2.2f | f = %2.4e | d(Z,X) %2.4e | rho % .4e\n', kiter, (frob(btdres(Y,U))/frob(Y))^2, dZ(end), rho_stable);

        mainloss_history = [mainloss_history; frob(btdres(Y,U))/frob(Y)];

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

        if (kiter > 5) && (mainloss_history(end) <= 1e-6) && (dZ(end) <= 1e-4)
            break;
        end
    end
    Y_hat = ful(U);
end
