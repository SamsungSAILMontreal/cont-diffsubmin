function [A, b, signal, noise] = gen_snr_problem(n_obs, signal_length, s, snr, stream)
    % Generate a random linear model b = A*signal + noise with A an
    % n_obs x signal_length matrix with iid Gaussian entries. Signal is of
    % length signal_length with entries in {-1,0,1} and sparsity s.
    % noise_level is measured in SNR_dB according to the standard formula 
    % (see https://sites.ualberta.ca/~msacchi/SNR_Def.pdf). 
    % For example, noise_level=20 will produce additive Gaussian
    % noise which has a signal-to-noise ratio (relative to measurements b)
    % of 20 dB. 
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


