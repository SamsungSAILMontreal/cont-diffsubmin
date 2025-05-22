function z_hat = babai(A, b, grid)
%%
%   compute the Babai estimation
%   find a sub-optimal solution for min_z ||R*z-y||_2
%   R - an upper triangular real matrix of n-by-n
%   y - a real vector of n-by-1
%   z_hat - resulting integer vector
%   A modified version of the code taken from
%   https://gist.github.com/sameerjagdale/216a9e47121ac5ca3dda
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

