clear all;clc;close all;
rng(2025)
addpath(genpath(pwd))
%% ------------------------------------------------------------------------
% This script is the main file for comparing LL1 algorithms
%%-------------------------------------------------------------------------
% Follows https://tensorlab.net/doc/ll1.html

%% Problem and tensor generation
size_tens = [10 11 12]*2;
L = [2 3 4]*2;
U = ll1_rnd(size_tens, L);
U
U{1}
U{1}{4}
% The CPD format can be requested using the 'OutputFormat' option:
options = struct;
options.OutputFormat = 'CPD';
U = ll1_rnd(size_tens, L, options);
% The result U is a cell of length R=3 with ∑rLr columns in the first two 
% factor matrices and R columns in the third factor matrix. 
% When using the CPD format, the parameter L is required for all algorithms, 
% as this parameter provides the necessary information on how the columns 
% are grouped.


% To expand the factorized representation U of the decomposition in multilinear 
% rank-(Lr,Lr,1) terms to a full tensor, ll1gen can be used:

size_tens = [10 11 12];
L = [2 3 4];
Ubtd = ll1_rnd(size_tens, L, 'OutputFormat', 'btd');
T    = ll1gen(Ubtd);
% T    = ll1gen(Ubtd, L); % L is optional for BTD format

% Ucpd = ll1_rnd(size_tens, L, 'OutputFormat', 'cpd');
% T    = ll1gen(U,L); % L is required for CPD format


%% Computing the decomposition in multilinear rank-(Lr,Lr,1) terms
% Our Solver
list_2fac = [];
nb_trials = 5;
for trial=1:nb_trials
    [T_hat, Uhat_2fac, mainloss_history, U0] = solver_2fac_ll1(T, L);
    list_2fac.Uhat{trial}=Uhat_2fac;
    list_2fac.lossfun{trial}=mainloss_history;
end


% Tensorlab
% [Uhat,output] = ll1(T, U0, L,'Display', 1, 'Initialization', init);
list_tensorlab = [];

for trial=1:nb_trials
    init = @ll1_rnd;
    [Uhat,output] = ll1(T, L,'Display', 1, 'Initialization', init);
    list_tensorlab.Uhat{trial}=Uhat;
    list_tensorlab.output{trial}=output;
end
%% 
% Find the best results among 5 runs of tensorlab and 2 fac updates
idx_best_tensorlab = 1;
idx_best_2fac = 1;
funval_best_tensorlab = list_tensorlab.output{1}.Algorithm.fval(end);
funval_best_2fac = list_2fac.lossfun{1}(end);
for trial=2:nb_trials
    if list_tensorlab.output{trial}.Algorithm.fval(end) < funval_best_tensorlab
        idx_best_tensorlab = trial;
        funval_best_tensorlab = list_tensorlab.output{trial}.Algorithm.fval(end);
    end

    if list_2fac.lossfun{trial}(end) < funval_best_2fac
        idx_best_2fac = trial;
        funval_best_2fac = list_2fac.lossfun{trial}(end);
    end

end

%% ------------------------------------------------------------------------
% Post-processing
%--------------------------------------------------------------------------
close all;
% idx_best_2fac = 2;
font_size = 15;
figure;
semilogy(1:list_tensorlab.output{idx_best_tensorlab}.Algorithm.iterations,sqrt(2*list_tensorlab.output{idx_best_tensorlab}.Algorithm.fval(2:end)),'-','LineWidth',2);
hold on
semilogy(1:length(list_2fac.lossfun{idx_best_2fac}),list_2fac.lossfun{idx_best_2fac}*frob(T),'-.','LineWidth',2);
text{1} = 'll1 - tensorlab';
text{2} = 'll1 - 2 Fac. Updates';
xlabel('iteration - $k$','Interpreter','latex','FontSize',font_size);
ylabel('$\| \mathcal{Y} - \sum_{r=1}^R \left(A_r B_r^T\right) \otimes c_r \|_F$',"Interpreter",'latex','FontSize',font_size);
legend(text,'Location','northwest','Orientation','horizontal',"Interpreter","latex",'FontSize',font_size)
grid on;

disp(['The relative Frobenius reconstruction with tensorlab error is ', num2str(frob(ful(list_tensorlab.Uhat{idx_best_tensorlab})-T)/frob(T))]);
disp(['The relative Frobenius reconstruction with ours error is ', num2str(frob(ful(list_2fac.Uhat{idx_best_2fac})-T)/frob(T))]);