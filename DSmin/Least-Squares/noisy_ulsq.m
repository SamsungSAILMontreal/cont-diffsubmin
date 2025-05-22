function [] = noisy_ulsq(n, m, use_gurobi, save_dir, task_id, job_id)
    % Recover random signal from {-1, 0, 2, 3}^(signal_length) using different
    % recovery approaches. Experiments are conduncted with varying snr levels
    % for the additive Gaussian noise.
    %
    % Inputs:
    %   n: The length of the random signal
    %
    %   m: Number of measurements
    %
    %   use_gurobi: boolean to indicate where gurobi should compute a
    %   global solution to the integer least squares problem.
    %
    %   save_dir: path to output file for writing results structure.
    %
    %   task_id: The task id corresponds to the number of measurements to use. For
    %            example:
    %            task_id = 1 -> SNR=15.
    %            task_id = 2 -> SNR=16.
    %            task_id = 3 -> SNR=17.
    %            .
    %            .
    %            .
    %            task_id = 16 -> SNR=30.
    %
    %   job_id: Job id number to record in output file.

    Folder = cd;
    Folder = fullfile(Folder, '..');
    addpath(Folder)
    
    mkdir(strcat(save_dir, "/", string(job_id)))
    setenv('MW_PCT_TRANSPORT_HEARTBEAT_INTERVAL', '6000')
    n_trials = 100;
    snrs = linspace(15, 30, 16);
    results = {};
    snr = snrs(task_id);
    parfor j = 1:n_trials
        grid = [-1, 0, 2, 3];
        k=length(grid);
        Problem = {};
        model = {};
        s = RandStream('Threefry');
        s.Substream = j;
        signal = grid(randi(s, 4, n, 1))';
        A = randn(s, m, n);
        unscaled_noise = randn(s, m, 1);
        
        y = A * signal;
        
        alpha = sqrt(10^(-snr/10) * (norm(y)/norm(unscaled_noise))^2);
        
        noise = alpha * unscaled_noise;
        
        b = y + noise;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Relax-and-Round solution %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tic;
        x_relax = lsqlin(A, b, zeros(n,n), zeros(n,1), zeros(n,n), zeros(n,1), min(grid) * ones(n,1), max(grid) * ones(n,1));
        x_rar = interp1(grid, grid, x_relax, 'nearest', 'extrap');
        time_rar(j) = toc;
        
        Q = A' * A;
        c = A' * b;
        
        Q_plus = max(Q, 0);
        Q_minus = Q - Q_plus;
        
        outer_iters = 50;
        inner_iters = 400;
        outer_eps = 1e-5;
        inner_eps = 1e-4;

        phi = @(x) grid(x + 1)';
        phi_inverse = @(x) sum((x == grid) .* [0:k-1], 2);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % DSmin solution %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tic;
        [xds_rar, ~, ds_info_rar] = dsmin_local_search(phi_inverse(x_rar), Q_plus, Q_minus, -c, grid, k, outer_iters, inner_iters, outer_eps, inner_eps);
        time = toc;
        xds_rar = phi(xds_rar);
        time_ds_rar(j) = sum(ds_info_rar.time);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Accelerated DSmin solution %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [xads_rar, ~, ads_info_rar] = dsmin_local_search(phi_inverse(x_rar), Q_plus, Q_minus, -c, grid, k, outer_iters, inner_iters, outer_eps, inner_eps, true);
        xads_rar = phi(xads_rar);
        time_ads_rar(j) = sum(ads_info_rar.time);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Old DSmin solution %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [xdsold_rar, ~, dsold_info_rar] = dsmin(phi_inverse(x_rar), Q_plus, Q_minus, -c, grid, k, outer_iters, inner_iters, outer_eps, inner_eps, []);
        xdsold_rar = phi(xdsold_rar);
        time_dsold_rar(j) = sum(dsold_info_rar.time);

        Problem.objective.P = Q;
        Problem.objective.q = -c;
        Problem.constraint.A = zeros(0, n);
        Problem.constraint.b = zeros(0, 1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % MIQP-ADMM solution %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tic;
        admm_history = miqp_admm(Problem, x_rar, grid, .5, outer_iters);
        time_admm(j) = toc;
        f_admm = inf;
        x_admm = zeros(n,1);
        for r = 1:outer_iters
            f_current = norm(A*admm_history(:, r) - b)^2;
            if f_current < f_admm
                f_admm = f_current;
                x_admm = admm_history(:, r);
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Gurobi (global) solution %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if use_gurobi && m >= n
            M = kron(sparse(eye(n)), [3, 1]);
            v = -ones(n, 1);
    
            encode = reshape([[0, 0]; [0, 1]; [1, 0]; [1, 1]], 1,4,2);
            binary_coding = tensorprod(double(signal==grid), encode, 2, 2);
            init_sol = reshape(permute(binary_coding, [3 1 2]), 2 * n, 1);
    
            model.Q = sparse(M' * Q * M);
            model.obj = 2*(v' * Q * M - c' * M);
            model.A = sparse(zeros(1, 2*n));
            model.rhs = 0;
            model.sense = '<';
            model.vtype = 'B';
            model.start = init_sol;
    
            sol = gurobi(model);
    
            x_gurobi = M * sol.x + v;
            time_gurobi(j) = sol.runtime;
            F_gurobi(j) = norm(A * x_gurobi - b)^2;
            correct_gurobi(j) = sum((x_gurobi - signal) == 0) / n;
            signals_gurobi(:, j) = x_gurobi;
        else
            F_signal(j) = norm(A * signal - b)^2;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Babai and OQP solution (if full-col rank)%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if m >= n
            tic;
            x_babai = babai(A, b, grid);
            time_babai(j) = toc;
            correct_babai(j) = sum((x_babai - signal) == 0) / n;
            F_babai(j) = norm(A * x_babai - b)^2;
            signals_babai(:, j) = x_babai;
        end
          
        tic;
        x_obq3 = obq(A, x_relax, grid, 2);
        time_obq3(j) = toc;
        
        correct_rar(j) = sum((x_rar - signal) == 0) / n;
        correct_admm(j) = sum((x_admm - signal) == 0) / n;
        correct_obq3(j) = sum(abs(x_obq3 - signal) <= 1e-4) / n;
        correct_ds_rar(j) = sum((xds_rar - signal) == 0) / n;
        correct_ads_rar(j) = sum((xads_rar - signal) == 0) / n;
        correct_dsold_rar(j) = sum((xdsold_rar - signal) == 0) / n;

        signals_rar(:, j) = x_rar;
        signals_admm(:, j) = x_admm;
        signals_obq3(:, j) = x_obq3;
        signals_ds_rar(:, j) = xds_rar;
        signals_ads_rar(:, j) = xads_rar;
        signals_dsold_rar(:, j) = xdsold_rar;
        
        F_rar(j) = norm(A * x_rar - b)^2;
        F_admm(j) = f_admm;
        F_obq3(j) = norm(A * x_obq3 - b)^2;
        F_ds_rar(j) = norm(A * xds_rar - b)^2;
        F_ads_rar(j) = norm(A * xads_rar - b)^2;
        F_dsold_rar(j) = norm(A * xdsold_rar - b)^2;
    end
    if use_gurobi && m >= n
        Fopt = F_gurobi;
    else
        Fopt = F_signal;
    end

    if m >= n
        results.correct_babai = correct_babai;
        results.optgap_babai = (F_babai - Fopt) ./ Fopt;
        results.time_babai = time_babai;
        results.signals_babai = signals_babai;
    end

    results.correct_rar = correct_rar;
    results.correct_admm = correct_admm;
    results.correct_obq3 = correct_obq3;
    results.correct_ds_rar = correct_ds_rar;
    results.correct_ads_rar = correct_ads_rar;
    results.correct_dsold_rar = correct_dsold_rar;

    results.signals_rar = signals_rar;
    results.signals_admm = signals_admm;
    results.signals_obq3 = signals_obq3;
    results.signals_ds_rar = signals_ds_rar;
    results.signals_ads_rar = signals_ads_rar;
    results.signals_dsold_rar = signals_dsold_rar;

    results.time_rar = time_rar;
    results.time_admm = time_admm;
    results.time_obq3 = time_obq3;
    results.time_ds_rar = time_ds_rar;
    results.time_ads_rar = time_ads_rar;
    results.time_dsold_rar = time_dsold_rar;

    results.optgap_rar = (F_rar - Fopt) ./ Fopt;
    results.optgap_admm = (F_admm - Fopt) ./ Fopt;
    results.optgap_obq3 = (F_obq3 - Fopt) ./ Fopt;
    results.optgap_ds_rar = (F_ds_rar - Fopt) ./ Fopt;
    results.optgap_ads_rar = (F_ads_rar - Fopt) ./ Fopt;
    results.optgap_dsold_rar = (F_dsold_rar - Fopt) ./ Fopt;

    if use_gurobi && m >= n
        results.correct_gurobi = correct_gurobi;
        results.time_gurobi = time_gurobi;
        results.signals_gurobi = signals_gurobi;
    end

    save(strcat(save_dir, "/", string(job_id), "/", string(snr)))
end