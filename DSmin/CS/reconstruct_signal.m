function [info] = reconstruct_signal(method, A, b, lambdas, true_signal, noise, time_all)
% RECONSTRUCT_SIGNAL  A function which attempts to reconstruct a sparse 
% signal using the technique referenced by method. A is the measurement
% matrix, true_signal is the sparse signal, and b is the set of
% measurements obtained by b = A * true_signal + noise. Lambdas is either a scalar
% or an array which changes the value of the sparsity penalty multiplier.
% This is used for methods which solve min_x ||Ax - b||^2 + lambda * F(x), 
% where supp(x) is the support of the vector x, and F is some sparsity 
% enforcing function . When lambdas is an array of values, then we loop over 
% the different values present and choose the reconstructed signal which most closely matches true_signal 
% in the L2 norm. For example, [best_signal, best_lsq_obj] = reconstruct_signal('lasso', A, b, [1,0.1,0.01], true_signal)
% attempts to reconstruct true_signal using lasso, looping over the different
% values of lambda = [1,0.1,0.01]. best_signal is the recovered signal
% which most closely matches true_signal, and best_lsq_obj is the least
% squares objective ||A * best_signal - b||^2.
    n = size(A, 2);
    fields = {'best_signal', 'best_lsq_obj', 'time', ...
                'f_min', 'f_vals', 'solutions', 'lambdas',... 
                'A', 'b', 'true_signal', 'best_supp_diff', 'noise',...
                'duality_gaps', 'local_min_gaps'};
    cells = cell(length(fields),1);
    info = cell2struct(cells,fields);
    %save copy of A, b, and true_signal
    info.A = A;
    info.b = b;
    info.noise = noise;
    info.true_signal = true_signal;
    Fvals = cell(1, length(lambdas));
    lambdas = sort(lambdas,'descend');
    temp = strsplit(method, '-');
    method_name = temp(1);
    warmstart = contains(method, "warmstart");
    lasso_init = contains(method, "LassoInit");
    if contains(method, "randomperm")
        break_ties = "random";
    else
        break_ties = "";
    end
    if strcmp(method_name, 'lasso')
        round_lattice = contains(method, "roundlattice");
        round_levelset = contains(method, "roundlevelset");
    end

    signal_vals = [-1, 0, 1];
    %Default "best" signal is the 0 signal
    best_signal = zeros(n, 1);
    best_lsq_obj = norm(b)^2;
    best_reconstruction_error = 1;
    best_supp_diff = sum(true_signal ~= 0);

    switch method_name
        
        case 'dsminLocalSearch'
            Q = A' * A;
            c = A' * b;
            % Lattice has values [0, 1, 2], so k = 3
            k = 3;
            Q_plus = max(Q, 0);
            Q_minus = Q - Q_plus;
%             F_quad = Quadratic(Q, -c, signal_vals);
%             G_quad = Quadratic(Q_minus, -c, signal_vals);
            %Use L0 norm, set q = 0
%             q = 0;
%             lq = Quasinorm(q, -1, 1, k);
%             H = Quadratic(-Q_plus, zeros(n,1), signal_vals);
            outer_iters = 25;
            inner_iters = 400;
            outer_eps = 1e-5;
            inner_eps = 1e-4;

            info.duality_gaps = cell(1, length(lambdas));
            info.local_min_gaps = [];
            info.times = cell(1, length(lambdas));

            signal2lattice = @(x) interp1(signal_vals, 0:(k-1), x, 'linear' , 'extrap');
            lattice2signal = @(x)interp1(0:(k-1), signal_vals, x, 'linear' , 'extrap');

            if lasso_init
                lasso_warmstart = @(x, lambda) signal2lattice(lasso_rounded(lattice2signal(x), A, b, lambda, signal_vals));
                reconstructor = @(x, lambda) dsmin_local_search(lasso_warmstart(x, lambda), Q_plus, Q_minus, -c, lambda, signal_vals, k, outer_iters, inner_iters, outer_eps, inner_eps, false, break_ties);
                init_sol = ones(n, 1);
                current_sol = init_sol;
            else
                %start from 0 vector (i.e. vector of all ones for the
                %lattice).
                reconstructor = @(x, lambda) dsmin_local_search(x, Q_plus, Q_minus, -c, lambda, signal_vals, k, outer_iters, inner_iters, outer_eps, inner_eps, false, break_ties);
                init_sol = ones(n, 1);
                current_sol = init_sol;
            end

        case 'dsmin'
            Q = A' * A;
            c = A' * b;
            % Lattice has values [0, 1, 2], so k = 3
            k = 3;
            Q_plus = max(Q, 0);
            Q_minus = Q - Q_plus;
%             F_quad = Quadratic(Q, -c, signal_vals);
%             G_quad = Quadratic(Q_minus, -c, signal_vals);
            %Use L0 norm, set q = 0
%             q = 0;
%             lq = Quasinorm(q, -1, 1, k);
%             H = Quadratic(-Q_plus, zeros(n,1), signal_vals);
            outer_iters = 25;
            inner_iters = 400;
            outer_eps = 1e-5;
            inner_eps = 1e-4;

            info.duality_gaps = cell(1, length(lambdas));
            info.local_min_check = cell(1, length(lambdas));
            info.local_min_gaps = [];
            info.times = cell(1, length(lambdas));

            signal2lattice = @(x) interp1(signal_vals, 0:(k-1), x, 'linear' , 'extrap');
            lattice2signal = @(x)interp1(0:(k-1), signal_vals, x, 'linear' , 'extrap');

            if lasso_init
                lasso_warmstart = @(x, lambda) signal2lattice(lasso_rounded(lattice2signal(x), A, b, lambda, signal_vals));
                reconstructor = @(x, lambda) dsmin(lasso_warmstart(x, lambda), Q_plus, Q_minus, -c, lambda, signal_vals, k, outer_iters, inner_iters, outer_eps, inner_eps, break_ties);
                init_sol = ones(n, 1);
                current_sol = init_sol;
            else
                %start from 0 vector (i.e. vector of all ones for the
                %lattice).
                reconstructor = @(x, lambda) dsmin(x, Q_plus, Q_minus, -c, lambda, signal_vals, k, outer_iters, inner_iters, outer_eps, inner_eps, break_ties);
                init_sol = ones(n, 1);
                current_sol = init_sol;
            end

        case 'lasso'
            % default maxiter in fista is 1000
            reconstructor = @(x, lambda) fista(@(y)0.5*norm(A*y-b,2)^2, @(y)A'*(A*y-b), @(y)norm(y,1), @(y,a)prox_l1(y,a), lambda, x);
            current_sol = zeros(n, 1);
            if round_lattice
                round_signal = @(x, lambda) interp1(signal_vals, signal_vals, x, 'nearest', 'extrap');
            elseif round_levelset
                round_signal = @(x, lambda) round_lasso_level_sets(x, A, b, lambda);
            else
                round_signal = @(x, lambda) x;
            end

        case 'omp'
            noise_level = norm(noise);
            sparsity = sum(true_signal ~= 0);
            extra_iters = round(sparsity/2);
            tic;
            [signals, ~] = OMP(A, b, sparsity, extra_iters, noise_level);
            info.time = toc;
    end

    if strcmp(method_name, 'omp') 
        n_signals = size(signals, 2);
        best_signal = zeros(n, 1);
        best_reconstruction_error = 1;
        best_supp_diff = sum(true_signal ~= 0);
        for i = 1:n_signals
            signal = interp1(signal_vals, signal_vals, signals(:, i), 'nearest', 'extrap');
            reconstruction_error = norm(signal - true_signal) / norm(true_signal);
            if reconstruction_error < best_reconstruction_error
                best_lsq_obj = norm(A*signal - b)^2;
                best_signal = signal;
                best_reconstruction_error = reconstruction_error;
            end
            supp_diff = norm((signal~=0) - (true_signal~=0), 1);
            if supp_diff < best_supp_diff
                best_supp_diff = supp_diff;
            end
        end

    else

        for i = 1:length(lambdas)
            lambda = lambdas(i);
            info.lambdas = [info.lambdas, lambda];
            if warmstart
                initialization = current_sol;
            else
                initialization = init_sol;
            end
            if strcmp(method_name, 'dsminLocalSearch') || strcmp(method_name, 'dsmin')
                [current_sol, f_min, ds_info] = reconstructor(initialization, lambda);
                Fvals{i} = ds_info.Fvals;
                info.time = [info.time, sum(ds_info.time)];
                info.times{i} = ds_info.time;
                info.duality_gaps{i} = ds_info.duality_gaps;
                info.local_min_gaps = [info.local_min_gaps, ds_info.local_min_gap];
                if strcmp(method_name, 'dsmin')
                    info.local_min_check{i} = ds_info.local_min_check;
                end
            else
                tic;
                [current_sol, f_min, f_vals] = reconstructor(initialization, lambda);
                info.time = [info.time, toc];
                Fvals{i} = f_vals;
            end
            info.solutions = [info.solutions, current_sol];
            info.f_min = [info.f_min, f_min];
            info.f_vals = Fvals;
            if strcmp(method_name, 'dsminLocalSearch') || strcmp(method_name, 'dsmin')
                signal = lattice2signal(current_sol);
                supp_diff = norm((signal~=0) - (true_signal~=0), 1);
            end
            if strcmp(method_name, 'pgm')
                signal = zeros(n ,1);
                lsq_sol = A(:, logical(current_sol)) \ b;
                signal(logical(current_sol)) = lsq_sol;
                supp_diff = norm((signal~=0) - (true_signal~=0), 1);
            end
            if strcmp(method_name, 'lasso')
                signal = round_signal(current_sol, lambda);
                supp_diff = norm((signal~=0) - (true_signal~=0), 1);
            end
            
            lsq_obj = norm(A*signal - b)^2;
            reconstruction_error = norm(signal - true_signal) / norm(true_signal);
    
            if reconstruction_error < best_reconstruction_error
                best_lsq_obj = lsq_obj;
                best_signal = signal;
                best_reconstruction_error = reconstruction_error;
            end

            if supp_diff < best_supp_diff
                best_supp_diff = supp_diff;
            end

            if reconstruction_error <= 10^-8 && ~time_all
                break; 
            end

        end
    end
    info.best_lsq_obj = best_lsq_obj;
    info.best_signal = best_signal;
    info.best_supp_diff = best_supp_diff;
end