function [xmin, Fmin, ds_info] = dsmin_local_search(x, Q_plus, Q_minus, c, lambda, grid, k, outer_iters, inner_iters, outer_eps, inner_eps, accelerate, break_ties)
% Approximates min_{x in {0, ..., k-1}^n} F(x) :=  G(x) - H(x)
% by applying DC algo (Algo 3 in the paper) to the equivalent continuous
% relaxation of min_{x in {0,...,k-1}^n} f(x) = 0.5*grid(x + 1)*Q*grid(x + 1)' + c^T*grid(x + 1)' + lambda * ||x||_0

% INPUT:
% x: Initial starting point x in {0,...,k-1}^n.
% Q_plus: Q_plus = max(Q,0).
% Q_minus: Q_minus = min(Q,0).
% grid: row vector of length (k-1) containing the (ordered) grid points.
% outer_iters: max number of outer iterations for DCA.
% inner_iters: max number of inner iterations for DCA.
% outer_eps: epsilon used in termination criteria for DCA.
% inner_eps: epsilon used in termination criteria for pairwise FW.
% break_ties: optional argument. to break ties when choosing permutations
% using a random permutation set break_ties="random".

% OUTPUT:
% xmin: The output from DCA.
% Fmin: The value of F at xmin.
% ds_info: struct containing additional information from each DCA
% iteration.

    if nargin < 12
        break_ties = [];
        accelerate = false;
    elseif nargin < 13
        break_ties = [];
    end

    q = 5;

    t = 1;

    n = size(x, 1);
    Q = Q_minus + Q_plus;
    F = @(x) 0.5 * grid(x + 1) * Q * grid(x + 1)' + grid(x + 1) * c + lambda*sum(grid(x+1) ~= 0);

    F_add = @(y, x, i) F_marginal(Q, c, lambda, y, x, grid, i, "add");
    F_rmv = @(y, x, i) F_marginal(Q, c, lambda, y, x, grid, i, "rmv");
    H_add = @(y, x, i) F_marginal(-Q_plus, zeros(n,1), 0, y, x, grid, i, "add");

    y0_G = 0.5 * grid(x + 1) * Q_minus * grid(x + 1)' + grid(x + 1) * c + lambda*sum(grid(x+1) ~= 0);
    y0_H = - 0.5 * grid(x + 1) * (Q_plus) * grid(x + 1)';
    y0_F = y0_G - y0_H;
 

    Fvals = zeros(outer_iters, 1); 
    F_current = F(x);
    Fvals(1) = F_current;
    x_current = x;
    duality_gaps = cell(1, outer_iters);
    time = zeros(outer_iters, 1);
    found_local_min = false;

    rho_current = theta_inverse(x_current, k, n);
    rho_prev = rho_current;

    ds_info_fields = {'Fvals', 'duality_gaps', 'time', 'local_min_gap'};
    ds_info_cells = cell(length(ds_info_fields),1);
    ds_info = cell2struct(ds_info_cells, ds_info_fields);
    ties = zeros(n, (k-1));
    for i=1:outer_iters
        tic;
        rho = rho_current;

        if strcmp(break_ties, "random")
            perm = reshape(randperm(n * (k-1)), n, k-1);
            ties = sort(perm, 2, "descend");
        end
        if accelerate
            t_next = (1 + sqrt(1 + 4 * t^2)) / 2;
            if i > 1
                z = rho_current + (t - 1) * (rho_current - rho_prev) / t_next;
                % [minz, maxz] = bounds(z, "all");
                if issorted(z, 2, 'descend')
                    [~, fz, ~, ~] = greedy_algorithm(z, y0_F, F_add, ties);
                        if fz <= max(Fvals(max(1, i-q):i))
                            rho = z;
                        end
                end
            end
            t = t_next;
        end

        [~, ~, diff] = get_best_neighbour(F_add, F_rmv, F_current, x_current, k);

        u = diff(1);
        if diff(2) > 0
            %we added a basis vector to x_best
            v = x_current(diff(1)) + 1;
            
            %This forces the new 1 at (u, v) to be ordered as the first 0 in rho.
            ties(u, v) = inf;
            W = greedy_algorithm(rho, y0_H, H_add, ties);

        elseif diff(2) < 0
            %we removed a basis vector to x_best
            v = x_current(diff(1));
            
            %This forces the new 0 at (u, v) to be ordered as the last 1 in rho
            ties(u, v) = -inf;
            W = greedy_algorithm(rho, y0_H, H_add, ties);

        else
            %we are already at a local min
            W = greedy_algorithm(rho, y0_H, H_add, ties);
            fprintf("DCA converged after %d iterations to a local min \n", i)
        end
        G_add = @(y, x, i) submod_major_marginal(Q_minus, c, lambda, -W, y, x, grid, i, "add");
%         mod_approx_H = ModFct(W);
%         submod_major = LatticeFctLinComb({G, mod_approx_H}, [1, -1]);
        [rho_new, dual_gaps_i] = fwpairwise(y0_G, G_add, rho, inner_iters, inner_eps);

        x_new = round_continuous_ext(rho_new, y0_G, G_add);
        duality_gaps{i} = dual_gaps_i;

        F_new = F(x_new);
        Fvals(i+1) = F_new;
        if F_current - F_new <= outer_eps
            time(i) = toc;
            found_local_min = true;
            % Fvals(i+1) = F_new;
            fprintf("DCA converged after %d iterations to a local min \n", i)
            break
        end
        time(i) = toc;
        rho_prev = rho_current;
        x_current = x_new;
        rho_current = theta_inverse(x_current, k, n);
        F_current = F_new;
        %Reset ties back to zero matrix if using default break ties
        ties(u, v) = 0;
    end

    ds_info.found_local_min = found_local_min;

    [~, F_best_neighbour] = get_best_neighbour(F_add, F_rmv, F_current, x_new, k, false);
    duality_gaps = duality_gaps(~cellfun(@isempty, duality_gaps));

    ds_info.local_min_gap = F_new - F_best_neighbour;
    ds_info.Fvals = Fvals(1:i+1);
    ds_info.duality_gaps = duality_gaps;
    ds_info.time = time(1:i);

    xmin = x_new;
    Fmin = F_new;
end