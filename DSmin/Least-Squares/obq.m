function [w] = obq(A, w, grid, variant)
    %An implementation of the Optimal Brain Quantization method given in:

    %Frantar, Elias, and Dan Alistarh. 
    %"Optimal brain compression: A framework for accurate post-training quantization and pruning." 
    % Advances in Neural Information Processing Systems 35 (2022): 4475-4488. 
   
    %Specifically, this code implements Algorithm 3 for obtaining a sub-optimal 
    %solution to min_{x in grid^n} ||A*x-b||_2
    %by quantizing a vector w to a point in grid^n.
    
    % INPUT:
    % A: mxn real matrix
    % grid: A 1d array of points in R^n.
    % w: an initial vector to quantize.

    % OUTPUT:
    % w: the quantized vector.
    [~,n] = size(A);
    H = 2 * (A' * A);
    if variant == 1
        H_inv = inv(H);
    else
        H_inv = inv(H + 0.01*mean(diag(H))*eye(n));
    end
    idxs = 1:n;
    for i = 1:n
        quant_weights = interp1(grid, grid, w(idxs), 'nearest', 'extrap');
        [~, j] = min((1 ./ diag(H_inv(idxs, idxs))) .* (quant_weights - w(idxs)).^2);
        p = idxs(j);
        q = quant_weights(j);
        w = w - H_inv(:, p) * (w(p) - q) ./ H_inv(p, p);
        %remove potential floating-point error by forcing w(p) to appropriate
        %grid value
        w(p) = q;
        H_inv = H_inv - H_inv(:, p) * H_inv(p, :) ./ H_inv(p,p);
        idxs = idxs(idxs ~= p);
    end
end

