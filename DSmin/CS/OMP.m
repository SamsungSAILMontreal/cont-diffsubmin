function [signals, residuals] = OMP(A, b, sparsity, extra_iters, noise_level)
    [m, n] = size(A);
    residuals = zeros(1, sparsity + extra_iters);
    signals = zeros(n, sparsity + extra_iters);
    switch nargin
        case 3
            extra_iters = 0;
            noise_level = 0;
        case 4
            noise_level = 0;
    end
    Q = zeros(m, n);
    %keep track of the support elements
    supp_indices = zeros(1, sparsity + extra_iters);
    %keep track of OMP iterations
    count = 0;
    r = b;
    while count < min(m, (sparsity + extra_iters)) && norm(r) >= noise_level
        %Break Conditions:
        %   1: The support has cardinality equal to the signal sparsity +   
        %      extra_iters.
        %   2: The support has cardinality equal to n_rows of A, hence
        %      we can already obtain a solution such that Ax = b.
        %   3: The magnitude of the residual ||r|| is less than 
        %      that magnitude of the the noise.

        %We proceed to update support set using OMP, following the approach
        %used by Stephen Becker *put reference here* 
        count = count + 1;
        [~, next_col] = max(abs(r' * A));
        supp_indices(count) = next_col;
        q_new = mgs_update(Q(:, 1:(count-1)), A(:, next_col));
        Q(:, count) = q_new;
        sol = Q(:, 1:count)' * b;
        r = b - Q(:, 1:count) * sol;

        %generate recovered signal
        signal = zeros(n, 1);
        x = A(:, supp_indices(1:count)) \ b;
        signal(supp_indices(1:count)) = x;
        signals(:, count) = signal;
        residuals(count) = norm(A * signal - b)^2;
    end
    signals = signals(:, 1:count);
    residuals = residuals(1:count);
end 