function history = miqp_admm(Problem, x_init, grid, rho, max_iter)
% A modified version of the code provided at https://github.com/cvxgrp/miqp_admm?tab=readme-ov-file which allows an initial point x_init to be given
% instead of using a random point.
% From the paper "A Simple Effective Heuristic for Embedded Mixed-Integer Quadratic Programming" by Reza takapoui, Nicholas Moehle, Stephen Boyd, and Alberto  %Bemporad, available at www.stanford.edu/~boyd/papers/miqp_admm.html
% Approximately solves the following problem using nonconvex ADMM
%
% minimize      (1/2) * x^T * P * x + q^T * x
% subject to    A*x = b
%               x(1: l1) == 0 or 1
%               x(l1 + 1: l1 + l2) == -1 or 0 or 1
%               x(l1 + l2  +1: l1 + l2 + l3) == -3 or -1 or 1 or 3
%               x(l1 + l2 + l3 + 1: l1 + l2 + l3 + l4) >= 0
%
% Inputs:
%   Problem:  A struct that contains objective and constraints,
%             objective is a struct containing P and q. constraints
%             is a struct containing A, b, l1, l2, l3, and l4.
%
%   rho:      The parameter rho is ADMM algorithm.
%
%   max_iter: Maximum number of iterations.
%
%
% Output:
%   history:  an n x max_iter matrix which contains the 
%             intermediate poins of the algorithm in its columns.  
%
% There is no guarantee that this solver finds the optimal solution.
% Tuning the input parameters can help the performance.

%%%%%%%%%%%%%%%%%%
% Default values %
%%%%%%%%%%%%%%%%%%
if nargin < 4
    repeat = 5;
end

if nargin < 3
    max_iter = 50;
end

if nargin < 2
    rho = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extracting problem information %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P = Problem.objective.P;
q = Problem.objective.q;
A = Problem.constraint.A;
b = Problem.constraint.b;
[m,n] = size(A);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Projection %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

projection = @(xx) interp1(grid, grid, xx, 'nearest', 'extrap');

%%%%%%%%%%%%%%%%%%%
% Preconditioning %
%%%%%%%%%%%%%%%%%%%
fprintf('Precondioning: ');
tic
q = q / norm(P);
P = P / norm(P);

E = 1 ./ norm(A',1)';
A = diag(E) * A;
b = E .* b;

b = b / norm(A);
A = A / norm(A);
Atb = A' * b;

M = inv(P + rho * (eye(n) + A' * A));
toc

%%%%%%%%%%%%%%%%%%%
% ADMM iterations %
%%%%%%%%%%%%%%%%%%%
rng('default')
history = zeros(n, max_iter);
tic

z = x_init;
u = zeros(n + m, 1);
for iter = 1:max_iter
    % x update
    x = M * (-q + rho * (z + Atb - A' * u(1:m)- u(m+1:end)));
    % z update
    z = projection(x + u(m+1:end));
    % u update
    u = u + [A * x; x] - [zeros(m,1); z]-[b; zeros(n,1)];
    % storing history
    history(:, iter) = z;
end
fprintf('.');

fprintf('\nSolving: ');
toc