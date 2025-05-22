%% Test add/remove
rng(1);
n = 20;
A = randn(100,n);
b =randn(100, 1);
lambda = 1;
grid = [-1, 0, 1];
f = @(x, lambda) 0.5*norm(A*grid(x+1)' - b)^2 - 0.5*norm(b)^2 + lambda*norm(grid(x+1))^2;
Q = A' * A + 2*lambda*eye(n);
c = -A' * b;
x = zeros(n,1);
y = f(x, lambda);
for i = 0:39
    j = mod(i,20) + 1;
    [y, x] = F_marginal(Q, c, y, x, grid, j, "add");
    y_true = f(x, lambda);
    if norm(y_true - y) >= 1e-10
        error("marginal gain function incorrect")
    end
end

for i = 0:39
    j = mod(i,20) + 1;
    [y, x] = F_marginal(Q, c, y, x, grid, j, "rmv");
    y_true = f(x, lambda);
    if norm(y_true - y) >= 1e-10
        error("marginal gain function incorrect")
    end
end

%% Test add/remove
rng(1);
n = 20;
A = randn(100,n);
b =randn(100, 1);
grid = [-1, 0, 2, 4];
f = @(x, lambda) 0.5*norm(A*grid(x+1)' - b)^2 - 0.5*norm(b)^2;
Q = A' * A;
c = -A' * b;
x = zeros(n,1);
y = f(x);
for i = 0:59
    j = mod(i,20) + 1;
    [y, x] = F_marginal(Q, c, y, x, grid, j, "add");
    y_true = f(x, lambda);
    if norm(y_true - y) >= 1e-10
        error("marginal gain function incorrect")
    end
end

for i = 0:59
    j = mod(i,20) + 1;
    [y, x] = F_marginal(Q, c, y, x, grid, j, "rmv");
    y_true = f(x, lambda);
    if norm(y_true - y) >= 1e-10
        error("marginal gain function incorrect")
    end
end
%% Test greedy algorithm implementation
rng(1);
n = 20;
A = randn(100, n);
b = randn(100, 1);
lambda = 1;
grid = [-1, 0, 1];
f = @(x, lambda) 0.5*norm(A*grid(x+1)' - b)^2 - 0.5*norm(b)^2 + lambda*norm(grid(x+1))^2;
Q = A' * A + 2*lambda*eye(n);
c = -A' * b;
x = zeros(n,1);
y = f(x, lambda);
rho = [0.5*rand(n,1) + 0.5, 0.5*rand(n,1)];
F_add = @(y, x, i) F_marginal(Q, c, y, x, grid, i, "add");
[w_true, f_true, Fmin_true] = bach_greedy_algorithm(rho, f, lambda);
[w, f, Fmin, ~] = greedy_algorithm(rho, y, F_add);
if norm(w_true - w) >= 1e-10
    error("greedy algorithm implementations don't match")
end

%% Test Pairwise Frank-Wolfe
rng(1);
n = 20;
maxiter = 100;
eps = 1e-5;
lambda = 1;
grid = [-1, 0, 1, 3];
Q = gen_random_submod_matrix(n);
f = @(x, params) 0.5*grid(x+1) * Q * grid(x+1)';
c = zeros(n,1);
x = zeros(n,1);
y = f(x, []);
rho = [0.5*rand(n,1) + 0.5, 0.5*rand(n,1)];
F_add = @(y, x, i) F_marginal(Q, c, y, x, grid, i, "add");
[rho_ours, gaps] = fwpairwise(y, F_add, rho, maxiter, eps);
[rho_true, gaps_true] = bach_fwpairwise(f, rho, [], maxiter, eps);
if norm(w_true - w) >= 1e-10
    error("greedy algorithm implementations don't match")
end