clear all; close all;clc;

addpath(genpath(pwd))

%% ------------------------------------------------------------------------
% Playground with tensorlab CPD basics
%%-------------------------------------------------------------------------

%%% Generation of random tensor

size_tens = [7 8 9]; R = 4;
% using a structure
options = struct;
% options.Real = @rand;
% options.Imag = @rand;
U = cpd_rnd(size_tens,R,options);

% using key-value pairs
U = cpd_rnd(size_tens, R, 'Real', @rand, 'Imag', @rand);

%% ------------------------------------------------------------------------
% Generate Tensors for CPD -  CHECK two-modes unfoldings
% -------------------------------------------------------------------------
disp('-------------------------------------------------------------------')
disp('---------------    CHECK two-modes unfoldings      ----------------')
disp('-------------------------------------------------------------------')
clear all; close all;clc;
R = 10;
szY = [15 16 17 18]; % You can uncomment the line below to use the original szY values
% szY = [8, 8, 8, 8];
% R = 4;
% szY = [10 15 10];

N = length(szY);
Factors = cell(1, N); % Preallocate a cell array for Factors

for n = 1:N
    Factors{n} = randn(szY(n), R); % Use randn for normally distributed random numbers
    % Factors{n} = rand(szY(n), R); % Uncomment this line to use uniformly distributed random numbers
end

weights_init = ones(R, 1); % Initialize weights as a column vector

% Y_hat
Y_hat.factors = Factors;
Y_hat.weights = ones(R,1);

% Y is a tensor of rank R
Y = cpdgen(Factors); 


%%% CHECK two-modes unfoldings 
modes = ndims(Y); % Get the number of modes (dimensions)
szY = size(Y); % Get the size of Y

for n = 1:(N-1)
    for m = (n+1):N
        
        modes_1 = sort([n, m]);
        modes_2 = setdiff(1:modes, modes_1); % Get the remaining modes

        % Yt = (Uk kr Ui)(Uj kr U1)^T != n, m
        sz2 = prod(szY(modes_2)); % Size of the unfolding for modes_2
        sz1 = prod(szY(modes_1)); % Size of the unfolding for modes_1

        % Unfolding mode-(n,m)    
        Y_nm = permute(Y, [modes_2, modes_1]); % Permute dimensions
        Y_nm = reshape(Y_nm, [sz2, sz1]); % Reshape to the new size

        A = kr(Y_hat.factors{modes_1(2)}, Y_hat.factors{modes_1(1)}); % Khatri-Rao product for modes_1
        Phi = kr(Y_hat.factors{modes_2(2)}, Y_hat.factors{modes_2(1)}); % Khatri-Rao product for modes_2

        % Calculate the norm difference
        norm_difference = norm(Y_nm - Phi * A'); % Using transpose for matrix multiplication
        disp(norm_difference); % Display the norm difference
    end
end

%% ------------------------------------------------------------------------
% Generate Tensors for CPD -  Check solver for the subproblem
% -------------------------------------------------------------------------
disp('-------------------------------------------------------------------')
disp('---------------   Check solver for the subproblem  ----------------')
disp('-------------------------------------------------------------------')
% clear all; close all;clc;
% load('matlab_matrix.mat')
% R = 10;
% szY = [15 16 17 18];
% N = length(szY);
% Factors = cell(1, N); % Preallocate a cell array for Factors
% Factors{1} = fac1; Factors{2} = fac2; Factors{3} = fac3; Factors{4} = fac4; 
% Y_hat.factors = Factors;
% Y_hat.weights = ones(R,1);
% Y_ref = Y;
% Y = cpdgen(Factors); 
% modes = ndims(Y);

%%% Check solver for the subproblem
n = 1; % MATLAB uses 1-based indexing
m = 2; % MATLAB uses 1-based indexing

modes_1 = sort([n, m]);
modes_2 = sort(setdiff(1:modes, modes_1)); % Get the remaining modes

sz2 = prod(szY(modes_2)); % Size of the unfolding for modes_2
sz1 = prod(szY(modes_1)); % Size of the unfolding for modes_1

% Unfolding mode-(n,m)    
Y_nm = permute(Y, [modes_2, modes_1]); % Permute dimensions
Y_nm = reshape(Y_nm, [sz2, sz1]); % Reshape to the new size

szXnm = [szY(n), szY(m), R]; % Size of X

% Initialize X0
U0 = Y_hat.factors{modes_1(1)}; % MATLAB uses 1-based indexing
V0 = Y_hat.factors{modes_1(2)}; % MATLAB uses 1-based indexing
X0 = kr(V0, U0 * diag(Y_hat.weights)); % Khatri-Rao product

% Prepare factors for the Khatri-Rao product
Factorx2 = cell(1, length(modes_2));
for i = 1:length(modes_2)
    Factorx2{i} = Y_hat.factors{modes_2(end - i + 1)}; % Reverse order
end
Phi = kr(Factorx2{:}); % Khatri-Rao product of the factors

rho = 1;
mu = 0;
maxiters = 100;

% Call the linreg_krp function (ensure this function is defined)
[X, Z, F1, F2, S, f, dZ, T, rhonew] = linreg_krp(Y_nm, Phi, rho, mu, szXnm, X0, maxiters, [], true);

% Plotting the results
figure;
loglog(f);
hold on;
loglog(dZ);
hold off;
xlabel('Iterations');
ylabel('Values');
title('Subsolver - Log-Log Plot of f and dZ ');
legend('f', 'dZ');
grid on;


%% ------------------------------------------------------------------------
% Generate Tensors for CPD -  Check main solver
% -------------------------------------------------------------------------
disp('-------------------------------------------------------------------')
disp('--------------- Check full solver - synthetic case ----------------')
disp('-------------------------------------------------------------------')
% Generate Y
N = length(szY);
Factors = cell(1, N); % Preallocate a cell array for Factors

for n = 1:N
    Factors{n} = randn(szY(n), R); % Use randn for normally distributed random numbers
    % Factors{n} = rand(szY(n), R); % Uncomment this line to use uniformly distributed random numbers
end

weights_init = ones(R, 1); % Initialize weights as a column vector

% Y_hat
Y_hat.factors = Factors;
Y_hat.weights = ones(R,1);

% Y is a tensor of rank R
Y = cpdgen(Factors); 


% Init for factors
Factors = cell(1, N); % Preallocate a cell array for Factors

for n = 1:N
    Factors{n} = randn(szY(n), R); % Use randn for normally distributed random numbers
    % Factors{n} = rand(szY(n), R); % Uncomment this line to use uniformly distributed random numbers
end

weights_init = ones(R, 1); % Initialize weights as a column vector

% Y_hat
Y_hat.factors = Factors;
Y_hat.weights = ones(R,1);

% Call of solver
maxoutiters = 100;
rho = 2;
maxiters = 20;
min_rho_stable = rho;
[Y_hat, mainloss_history] = solver_2fac_CPD(Y,R,Y_hat,rho,mu,maxoutiters,maxiters,min_rho_stable);

% Computation of cp_sensitivity
Y_hat.shape = szY;
Y_hat.rank = R;
disp(['The cp_sensitivity is ', num2str(cp_sensitivity(Y_hat))]);

% post-processing (tensorlab does not deal with weights)
Y_hat.factors{1} = Y_hat.factors{1} * diag(Y_hat.weights);
Y_hat.weights = ones(R,1);
Y_hat_full = cpdgen(Y_hat.factors); 
disp(['The relative Frobenius reconstruction error is ', num2str(frob(Y - Y_hat_full)/frob(Y))]);

% Plotting the results
figure;
loglog(mainloss_history);
hold on;
xlabel('Iterations');
ylabel('Values');
title('Full solver - Log-Log Plot of err');
legend('f');
grid on;


