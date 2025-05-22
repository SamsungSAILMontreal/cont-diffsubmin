function [xmin, Fmin, ds_info] = dsmin(x, Q_plus, Q_minus, c, grid, k, outer_iters, inner_iters, outer_eps, inner_eps, break_ties)
% Approximates min_{x in {0, ..., k-1}^n} F(x) :=  G(x) - H(x)
% by applying DC algo to the equivalent continuous relaxation problem
% min_{rho in [0,1]↓^{n x (k-1)}} f(rho) := g(rho) - h(rho)
% INPUT:
% OUTPUT:
    n = size(x, 1);
    Q = Q_minus + Q_plus;
    F = @(x) 0.5 * grid(x + 1) * Q * grid(x + 1)' + grid(x + 1) * c;

    F_add = @(y, x, i) F_marginal(Q, c, y, x, grid, i, "add");
    F_rmv = @(y, x, i) F_marginal(Q, c, y, x, grid, i, "rmv");
    H_add = @(y, x, i) F_marginal(-Q_plus, zeros(n,1), y, x, grid, i, "add");

    y0_G = 0.5 * grid(x + 1) * Q_minus * grid(x + 1)' + grid(x + 1) * c;
    y0_H = - 0.5 * grid(x + 1) * (Q_plus) * grid(x + 1)';

    F_current = F(x);
    Fvals = zeros(outer_iters + 1, 1);
    % Fvals(1) = F_current;
    duality_gaps = cell(1, outer_iters);
    time = zeros(outer_iters, 1);
    converged = false;
    local_min_check = zeros(outer_iters, 1);

    ds_info_fields = {'Fvals', 'duality_gaps', 'time', 'local_min_gap', 'local_min_check'};
    ds_info_cells = cell(length(ds_info_fields),1);
    ds_info = cell2struct(ds_info_cells, ds_info_fields);

    if isnan(x(1))
        rho = fliplr(cumsum(rand(n,k-1),2));
        rho = rho - min(rho(:));
        rho = rho / max(rho(:));
        [~, F_old] = greedy_algorithm(rho, F);
    elseif isinf(x(1))
        x = randi([0, k-1], n, 1);
        rho = theta_inverse(x, k, n);
        F_old = F(x);
    else
        rho = theta_inverse(x, k, n); 
        F_old = F_current;
    end
    Fvals(1) = F_old;
    for i=1:outer_iters
        tic;
        if strcmp(break_ties, "random")
            perm = reshape(randperm(n * (k-1)), n, k-1);
            ties = sort(perm, 2, "descend");
        else
            ties = [];
        end
        W = greedy_algorithm(rho, y0_H, H_add, ties);
        G_add = @(y, x, i) submod_major_marginal(Q_minus, c, -W, y, x, grid, i, "add");
        [rho, dual_gaps_i] = fwpairwise(y0_G, G_add, rho, inner_iters, inner_eps);

        x_new = round_continuous_ext(rho, y0_G, G_add);
        duality_gaps{i} = dual_gaps_i;

        F_new = F(x_new);

        if F_old - F_new <= outer_eps
            local_min_check(i) = 1;
            [x_best_neighbour, F_best_neighbour] = get_best_neighbour(F_add, F_rmv, F_new, x_new, k, false);
            if F_best_neighbour < F_new
                x_new = x_best_neighbour;
                F_new = F_best_neighbour;
                rho = theta_inverse(x_new, k, n);
            else
                ds_info.local_min_gap = F_new - F_best_neighbour;
                time(i) = toc;
                converged = true;
                Fvals(i+1) = F_new;
                fprintf("DCA converged after %d iterations to a local min \n", i)
                break
            end
        else
            local_min_check(i) = 0;
        end
        time(i) = toc;
        Fvals(i+1) = F_new;
        F_old = F_new;
    end

    ds_info.converged = converged;
    if ~converged
        [~, F_best_neighbour] = get_best_neighbour(F_add, F_rmv, F_new, x_new, k, false);
    end
    duality_gaps = duality_gaps(~cellfun(@isempty, duality_gaps));

    ds_info.Fvals = Fvals(1:(i+1));
    ds_info.local_min_gap = F_new - F_best_neighbour;
    ds_info.duality_gaps = duality_gaps;
    ds_info.time = time(1:i);
    ds_info.local_min_check = local_min_check(1:i);

    xmin = x_new;
    Fmin = F_new;
end