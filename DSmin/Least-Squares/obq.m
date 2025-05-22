function [w] = obq(A, w, grid, variant)
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

