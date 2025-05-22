function [A, b, signal, noise] = gen_snr_problem(n_obs, signal_length, s, snr, stream)
    %generate the random signal first for reproducibility
    supp = randperm(stream, signal_length, s);
    r = (rand(stream, 1, s)<.5)*2 - 1;
    signal = zeros(signal_length, 1);
    signal(supp) = r;

    A = randn(stream, n_obs, signal_length) ./ sqrt(n_obs);
    unscaled_noise = randn(stream, n_obs, 1);
    y = A * signal;

    alpha = sqrt(10^(-snr/10) * (norm(y)/norm(unscaled_noise))^2);

    noise = alpha * unscaled_noise;

    b = y + noise;
end


