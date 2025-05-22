function [A, b, signal, noise] = gen_random_signal_problem(n_obs, signal_length, s, noise_level, seed)
    %generate the random signal first for reproducibility
    supp = randperm(seed, signal_length, s);
    r = (rand(seed, 1, s)<.5)*2 - 1;
    signal = zeros(signal_length, 1);
    signal(supp) = r;

    rand_mat = randn(seed, n_obs, signal_length + 1);
    A = 1/(sqrt(n_obs))*rand_mat(:, 1:signal_length);
    unscaled_noise = rand_mat(:, signal_length+1);

    noise = noise_level * unscaled_noise;
    b = A * signal + noise;
end

