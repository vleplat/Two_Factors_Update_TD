clear; clc; close all;

% Focused validation for the two BTD examples reported in the paper.
% It compares only the ADMM block order, using the same tensor and the same
% initialization in each paired run:
%   revised: X -> Z -> T
%   reference: Z -> X -> T
%
% Tensorlab is optional because the immediate purpose is to measure the
% impact of the order exchange. Set run_tensorlab=true to include it.

rng(2025, 'twister');
addpath(genpath(pwd));

nb_trials = 5;
run_tensorlab = false;
mu = 0;
maxoutiters = 200;
maxiters = 100;
min_rho_stable = 0.1;

cases = struct([]);
cases(1).name = 'small';
cases(1).size_tens = [10 11 12];
cases(1).L = [2 2 3];
cases(1).rho = 2;

cases(2).name = 'mid';
cases(2).size_tens = [50 55 60];
cases(2).L = [6 9 12];
cases(2).rho = 10;

results_root = fullfile(pwd, 'Results', 'Revision_Order_Check');
if ~exist(results_root, 'dir')
    mkdir(results_root);
end

for icase = 1:numel(cases)
    cfg = cases(icase);
    fprintf('\n============================================================\n');
    fprintf('BTD case %s: size [%s], L = [%s]\n', cfg.name, ...
        num2str(cfg.size_tens), num2str(cfg.L));
    fprintf('============================================================\n');

    % Generate one exact input tensor for this case.
    Utrue = ll1_rnd(cfg.size_tens, cfg.L, 'OutputFormat', 'btd');
    Tdata = ll1gen(Utrue);

    results = struct();
    results.config = cfg;
    results.Utrue = Utrue;

    for trial = 1:nb_trials
        fprintf('\n--- trial %d/%d ---\n', trial, nb_trials);
        Uinit = ll1_rnd(cfg.size_tens, cfg.L, 'OutputFormat', 'btd');
        results.Uinit{trial} = Uinit;

        [~, U_xfirst, hist_xfirst] = solver_2fac_ll1( ...
            Tdata, cfg.L, Uinit, cfg.rho, mu, maxoutiters, maxiters, ...
            min_rho_stable, 'RAND');
        results.xfirst.U{trial} = U_xfirst;
        results.xfirst.history{trial} = hist_xfirst;
        results.xfirst.final_relerr(trial) = ...
            frob(btdres(Tdata, U_xfirst)) / frob(Tdata);

        [~, U_zfirst, hist_zfirst] = solver_2fac_ll1_zfirst_reference( ...
            Tdata, cfg.L, Uinit, cfg.rho, mu, maxoutiters, maxiters, ...
            min_rho_stable, 'RAND');
        results.zfirst.U{trial} = U_zfirst;
        results.zfirst.history{trial} = hist_zfirst;
        results.zfirst.final_relerr(trial) = ...
            frob(btdres(Tdata, U_zfirst)) / frob(Tdata);

        if run_tensorlab
            [U_tlab, output_tlab] = ll1(Tdata, Uinit, 'Display', 1);
            results.tensorlab.U{trial} = U_tlab;
            results.tensorlab.output{trial} = output_tlab;
            results.tensorlab.final_relerr(trial) = ...
                frob(btdres(Tdata, U_tlab)) / frob(Tdata);
        end
    end

    [~, best_x] = min(results.xfirst.final_relerr);
    [~, best_z] = min(results.zfirst.final_relerr);

    fprintf('\nFinal relative errors (X-first):\n');
    disp(results.xfirst.final_relerr);
    fprintf('Final relative errors (Z-first reference):\n');
    disp(results.zfirst.final_relerr);

    figure;
    semilogy(results.xfirst.history{best_x}, '-.', 'LineWidth', 2);
    hold on;
    semilogy(results.zfirst.history{best_z}, '--', 'LineWidth', 2);
    grid on;
    xlabel('outer iteration k');
    ylabel('relative Frobenius error');
    legend('X first, then Z', 'Z first, then X', 'Location', 'best');
    title(sprintf('BTD order comparison: %s case', cfg.name));

    case_dir = fullfile(results_root, cfg.name);
    if ~exist(case_dir, 'dir')
        mkdir(case_dir);
    end
    save(fullfile(case_dir, 'order_comparison.mat'), ...
        'results', 'Tdata', '-v7.3');
    savefig(fullfile(case_dir, 'order_comparison.fig'));
    print(gcf, fullfile(case_dir, 'order_comparison.eps'), '-depsc');
end
