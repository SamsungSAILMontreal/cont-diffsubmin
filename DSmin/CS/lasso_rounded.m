function [sol_rounded] = lasso_rounded(x_init, A, b, lambda, round_values)
    %Solves the problem 
    %min_{x in R^n} ||Ax-b||^2 + lambda * ||x||_1 (P)
    %using the fista solver from https://www.tau.ac.il/~becka/solvers/fista.
    %Each element in the solution x_min of (P) is then rounded to the nearest value in 
    % round_values.
    sol = fista(@(y)0.5*norm(A*y-b,2)^2, @(y)A'*(A*y-b), @(y)norm(y,1), @(y,a)prox_l1_box(y,a), lambda, x_init);
    sol_rounded = interp1(round_values, round_values, sol, 'nearest', 'extrap');
end