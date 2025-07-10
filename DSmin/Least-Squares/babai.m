function z_hat = babai(A, b, grid)
%   A modified version of the code taken from
%   https://gist.github.com/sameerjagdale/216a9e47121ac5ca3dda.

%   compute the box-constrained Babai point to obtain a sub-optimal 
%   solution for min_{x in grid^n} ||A*x-b||_2
%   see, for example, https://www.cs.mcgill.ca/~chang/pub/WenCT17.pdf
%   Section II for more details.
% INPUT:
% A: mxn real matrix
% b: m-dimensional array of measurements
% grid: A 1d array of points in R^n.

% OUTPUT:
% z_hat: the box-constrained Babai point.
    [Q, R] = qr(A, "econ");
    y = Q' * b;
    n = length(y);
    z_hat = zeros(n,1);
    z_hat(n) = interp1(grid, grid, y(n)./R(n, n), 'nearest', 'extrap');
    for k=n-1:-1:1
        par=R(k,k+1:n)*z_hat(k+1:n);
        ck=(y(k)-par)./R(k,k);
        z_hat(k)=interp1(grid, grid, ck, 'nearest', 'extrap');
    end
end

