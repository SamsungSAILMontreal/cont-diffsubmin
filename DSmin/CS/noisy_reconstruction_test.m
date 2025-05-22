function [results] = noisy_reconstruction_test(signal_length, s, job_id, task_id, resultsdir)
    % Recover random s-sparse signal from {-1,0,1}^(signal_length) using different
    % recovery approaches. Experiments are conduncted over a range of
    % snr levels from 0dB up to 20dB.
    %
    % Inputs:
    %   signal_length: The length of the random signal
    %
    %   s: Sparsity level (i.e. number of nonzero elements).
    %
    %   job_id: Job id number to record in output file.
    %
    %   task_id: The task id corresponds to the method you wish to run. For
    %            example:
    %            task_id = 1 -> dsminLocalSearch-warmstart.
    %            task_id = 2 -> dsminLocalSearch-warmstart-LassoInit.
    %            task_id = 3 -> dsminLocalSearch-randomperm-warmstart.
    %            task_id = 4 -> dsminLocalSearch-randomperm-warmstart-LassoInit.
    %            task_id = 5 -> lasso-roundlattice-warmstart.
    %            task_id = 6 -> omp.
    %
    %   resultsdir: path to output file for writing results structure. 
    %
    % Output:
    %   results:  A structure containing the recovered signals and corresponding
    %             metrics.
    
    disp([class(s), class(job_id)]);

    Folder = cd;
    Folder = fullfile(Folder, '..');
    addpath(Folder)

    methods = [ "dsminLocalSearch-warmstart",...
                "dsminLocalSearch-warmstart-LassoInit",...
                "dsminLocalSearch-randomperm-warmstart",...
                "dsminLocalSearch-randomperm-warmstart-LassoInit",...
                "lasso-roundlattice-warmstart",...
                "omp",...
                "dsmin-warmstart",...
                "dsmin-warmstart-LassoInit"];

    method = methods(task_id);

    %experiment hyperparameters
    n_trials = 100;
    n_obs = 128;
    snr_levels = linspace(0, 20, 21);

    addpath(genpath(pwd));
    [~, git_hash_string] = system('git rev-parse --short HEAD');
    dirname = sprintf(resultsdir + "/signal_length=%d,sparsity=%d,n_trials=%d,%s-%s", signal_length, s, n_trials, erase(git_hash_string,newline), num2str(job_id));
%     dirname = sprintf(resultsdir + "/signal_length=%d,sparsity=%d,n_trials=%d,jobid=%s", signal_length, s, n_trials, num2str(job_id));
    mkdir(dirname);

    results = {};
    
    results.method = method;
    results.n_obs = n_obs;
    results.sparsity = s;
    results.signal_length = signal_length;
    results.n_trials = n_trials;

    lambdas = [1, 0.1, 0.01, 0.001, 0.0001, 0.00001];

    for i = 1:length(snr_levels)
        snr = snr_levels(i);
        l2error = zeros(1, n_trials);
        lsq_objective = zeros(1, n_trials);
        sparsity = zeros(1, n_trials);
        suppdiff = zeros(1, n_trials);
        opt_info = cell(1, n_trials);
        parfor j = 1:n_trials
            stream = RandStream('Threefry');
            stream.Substream = j;
            [A, b, signal, noise] =  gen_snr_problem(n_obs, signal_length, s, snr, stream);
            info = reconstruct_signal(method, A, b, lambdas, signal, noise, false);
            l2error(j) = norm(info.best_signal - signal) / norm(signal);
            lsq_objective(j) = info.best_lsq_obj;
            sparsity(j) = sum(info.best_signal ~= 0);
            suppdiff(j) = info.best_supp_diff;
            opt_info{j} = info;
        end
        results.((strcat('snr',num2str(snr)))).l2error = l2error;
        results.((strcat('snr',num2str(snr)))).lsq_objective = lsq_objective;
        results.((strcat('snr',num2str(snr)))).sparsity = sparsity;
        results.((strcat('snr',num2str(snr)))).suppdiff = suppdiff;
        results.((strcat('snr',num2str(snr)))).opt_info = opt_info;
    end

    %% Save results
    filename = dirname + sprintf("/%s", datetime('now'));
    fprintf("saving results to %s", filename)
    save(filename);
end