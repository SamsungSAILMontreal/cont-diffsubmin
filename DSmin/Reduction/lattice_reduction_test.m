% Minimize 0.5 x^T Q x where q_ij <= 0 for i != j and x in the
% discretized hypercube [lb, ub]^n. Here, we compare minimization
% using a lattice-submodular function and the corresponding
% set function reduction given in Bach, Francis. "Submodular functions: from discrete to continuous domains." 
% Mathematical Programming 175 (2019), Section 4.4.


% where to save results
dirname = '/home/georgeao/scratch/ReductionResults';
% mkdir(dirname);

%ambient dimension n
n = 50;

%lower and upper bounds for the n-dimensionsal cube
lb = -1;
ub = 1;

param.lb = lb;
param.n = n;

%define number of Frank-Wolfe iterations
n_iters = 200;

%this prevents the Frank-Wolfe algorithm from terminating early if the 
%duality gap is small
eps = 0;

%define different discretization levels k, together with the number of
%trials to perform for each discretization level
k_vals = [10, 20, 30, 40, 50];
n_trials = 50;

seeds = 1:n_trials;

results = {};
fields = {'avg_gap', 'avg_time'};
cells = cell(length(fields),1);

for k = 1:length(k_vals)
    k_val = num2str(k_vals(k));
    results.((strcat('k_val',k_val))).('lattice') = cell2struct(cells,fields);
    results.((strcat('k_val',k_val))).('reduction') = cell2struct(cells,fields);
    results.((strcat('k_val',k_val))).('lattice').('dual_gaps') = cell(1, n_trials);
    results.((strcat('k_val',k_val))).('reduction').('dual_gaps') = cell(1, n_trials);
    results.((strcat('k_val',k_val))).('lattice').('relative_objective_gaps') = zeros(1, n_trials);
    results.((strcat('k_val',k_val))).('reduction').('relative_objective_gaps') = zeros(1, n_trials);
end

for i = 1:length(k_vals)
    k = k_vals(i);
    param.k = k;

    %start both approaches at the same point i.e. the zero matrix
    zero_matrix_lattice = zeros(n, k-1);
    zero_matrix_reduction = zeros(n * (k-1), 1);

    F_reduction = @(x, param) lattice_reduction_quadratic(x, param);
    F_lattice = @(x, param) lattice_quadratic(x, param);

    param.m = (ub-lb)/(k-1);

    avg_gap_lattice = 0;
    avg_gap_reduction = 0;

    times_lattice = zeros(1, n_trials);
    times_reduction = zeros(1, n_trials);

    for j = 1:n_trials
        rng(seeds(j))
        param.Q = gen_random_submod_matrix(n);
        param.L = max(sum(abs(param.Q), 2));
    
        %Time the lattice minimization problem
        tic;
        [rho_lattice, lattice_gaps] = bach_fwpairwise(F_lattice, zero_matrix_lattice, param, n_iters, eps);
        times_lattice(j) = toc;
        % avg_time_lattice = avg_time_lattice + time_lattice;
        x_lattice = bach_theta_minimizer(rho_lattice, F_lattice, param);
        avg_gap_lattice = avg_gap_lattice + min(lattice_gaps);
        results.((strcat('k_val',num2str(k)))).('lattice').('dual_gaps'){j} = lattice_gaps;
        
        %Time the reduction minimization problem
        tic;
        [rho_reduction, reduction_gaps] = bach_fwpairwise(F_reduction, zero_matrix_reduction, param, n_iters, eps);
        times_reduction(j) = toc;
        % avg_time_reduction = avg_time_reduction + toc;
        rho_reduction = reshape(bach_theta_minimizer(rho_reduction, F_reduction, param), n, k-1);
        x_reduction = sum(rho_reduction, 2);
        avg_gap_reduction = avg_gap_reduction + min(reduction_gaps);
        results.((strcat('k_val',num2str(k)))).('reduction').('dual_gaps'){j} = reduction_gaps;
        
        %To ensure that the solution rho_reduction is in {0,1}↓^{n x (k-1)}
        %use map_row_noninc and check that output is identical
        rho_mapped = map_row_noninc(rho_reduction);
        if isequal(rho_mapped, rho_reduction)
            %We can map rho_reduction to a vector x_reduction in the
            %lattice
            lattice_opt_val = F_lattice(x_lattice, param);
            reduction_opt_val = F_lattice(x_reduction, param);
            best_val = min(lattice_opt_val, reduction_opt_val);
            results.((strcat('k_val',num2str(k)))).('lattice').('relative_objective_gaps')(j) = abs((lattice_opt_val - best_val) / best_val);
            results.((strcat('k_val',num2str(k)))).('reduction').('relative_objective_gaps')(j) = abs((reduction_opt_val - best_val) / best_val);
        else
            lattice_opt_val = F_lattice(x_lattice, param);
            results.((strcat('k_val',num2str(k)))).('lattice').('relative_objective_gaps')(j) = 0;
            results.((strcat('k_val',num2str(k)))).('reduction').('relative_objective_gaps')(j) = nan;
        end
    end
    avg_time_lattice = mean(times_lattice);
    avg_time_reduction = mean(times_reduction);
    avg_gap_lattice = avg_gap_lattice / n_trials;
    avg_gap_reduction = avg_gap_reduction / n_trials;

    results.((strcat('k_val',num2str(k)))).('lattice').('avg_gap') = avg_gap_lattice;
    results.((strcat('k_val',num2str(k)))).('reduction').('avg_gap') = avg_gap_reduction;

    results.((strcat('k_val',num2str(k)))).('lattice').('avg_time') = avg_time_lattice;
    results.((strcat('k_val',num2str(k)))).('reduction').('avg_time') = avg_time_reduction;
end

filename = dirname + sprintf("/%d-%s", seeds(1), datetime('now'));
fprintf("saving results to %s", filename)
save(filename);

%% plot graphics
% results_dir = '/Users/georgeorfanides/Documents/GitHub/diff-sub-min/Code/DS-Lattice/results/ReductionTest/';
% file_name = '1-01-Nov-2024 19:37:59.mat';
% load(strcat(results_dir, file_name))
% 
% 
% 
% fields = string(fieldnames(results));
% k_vals = [];
% 
% for i = 1:length(fields)
%     k_vals(i) = str2num(erase(fields(i), 'k_val'));
% end
% 
% length_k = length(k_vals);
% time_vals = zeros(length_k, 2);
% dual_vals = zeros(length_k, 2);
% dual_std = zeros(length_k, 2);
% relative_vals = zeros(length_k, 1);
% relative_std = zeros(length_k, 1);
% 
% for i = 1:length_k
%     time_vals(i, 1) = results.(fields(i)).lattice.avg_time;
%     time_vals(i, 2) = results.(fields(i)).reduction.avg_time;
% 
%     dual_vals(i, 1) = results.(fields(i)).lattice.avg_gap;
%     dual_vals(i, 2) = results.(fields(i)).reduction.avg_gap;
% 
%     relative_vals(i) = mean(results.(fields(i)).reduction.relative_objective_gaps);
%     relative_std(i) = std(results.(fields(i)).reduction.relative_objective_gaps);
% end
% 
% clf;
% plot([k_vals', k_vals'] , time_vals, 'linewidth', 2, 'marker', '*', 'markersize',10);
% set(gca,'fontsize',18)
% set(gca, 'YScale', 'log')
% legend('Lattice', 'Reduction', 'Location','northeastoutside')
% title('Average Time vs Discrectization Level k')
% xlabel('k') 
% ylabel('Average Time')
% saveas(gcf, strcat(results_dir,'avg_time.fig'))
% exportgraphics(gcf,strcat(results_dir,'avg_time.pdf'))
% 
% clf;
% errorbar(k_vals', relative_vals, relative_std, 'linewidth', 2, 'marker', '*', 'markersize',10);
% set(gca,'fontsize',18)
% title('Average Relative Gap vs Discrectization Level k')
% xlabel('k') 
% ylabel('Average Relative Gap')
% saveas(gcf, strcat(results_dir,'avg_relative_gap.fig'))
% exportgraphics(gcf,strcat(results_dir,'avg_relative_gap.pdf'))
% 
% clf;
% plot([k_vals', k_vals'] , dual_vals, 'linewidth', 2, 'marker', '*', 'markersize',10);
% set(gca,'fontsize',18)
% set(gca, 'YScale', 'log')
% legend('Lattice', 'Reduction', 'Location','northeastoutside')
% title('Average Duality Gap vs Discrectization Level k')
% xlabel('k') 
% ylabel('Average Duality Gap')
% saveas(gcf, strcat(results_dir,'avg_dual_gap.fig'))
% exportgraphics(gcf,strcat(results_dir,'avg_dual_gap.pdf'))