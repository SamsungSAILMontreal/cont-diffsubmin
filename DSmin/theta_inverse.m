function [rho] = theta_inverse(x, k, n)
    % map x in [0:k-1]^n to a row-nonincreasing binary matrix rho in
    % {0,1}↓^{n x (k-1)}
    rho = (1: k-1) .* ones(n,1) <= x;
end

