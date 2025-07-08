function [w, f, Fmin, xmin] = greedy_algorithm(rho, y0, F_add, ties)
%A modified version of the greedy algorithm from:

% Bach, Francis. "Submodular functions: from discrete to continuous domains." 
% Mathematical Programming 175 (2019): 419-459. 
% (http://www.di.ens.fr/~fbach/submodular_multi_online.zip) Retrieved January 22, 2024.

% For submodular F, greedy_algorithm(rho, y0, F_add, ties) computes a
% subgradient of the continuous extension of F at rho.

% INPUT:
% rho: An n*(k-1) row-nonincreasing matrix with elements in [0,1].
% y0: y0 = F(0).
% F_add: Given initial y = F(x), F_add computes y_add = F(x + e_i).
% ties: Optional argument. to break ties when choosing permutations
% using a random permutation set break_ties="random".

% OUTPUT:
% w: An n*(k-1) row-nonincreasing matrix which is the subgradient of the 
% continuous extension of F at rho.
% f: The value of the continuous extension at rho. 
% xmin: Round rho to a point xmin in {0,...,k-1}^n such that F(xmin) is less than
% or equal to the continuous extension of F at rho.
% Fmin: F(xmin).

[n, k] = size(rho);
w = zeros(size(rho));
k = k + 1;

% first order all rhos (does preserve the ordering within rows if equal
% values)
if nargin < 4 || isempty(ties)
    [~, s] = sort( reshape(rho',n*(k-1),1), 'descend');
    [js, is] = ind2sub([k-1 n],s);
else
    [~, s] = sortrows([reshape(rho',n*(k-1),1), reshape(ties',n*(k-1),1)], 'descend');
    [js, is] = ind2sub([k-1 n],s);
end

% now go through all elements
xold = zeros(n,1);
Fold = y0;
H0 = Fold;

xmin = xold;
Fmin = Fold;

for i=1:n*(k-1)
    [Fnew, xnew] = F_add(Fold, xold, is(i));
    if (Fnew<Fmin), xmin = xnew; Fmin = Fnew; end
    w(is(i),js(i)) = Fnew - Fold;
    xold = xnew;
    Fold = Fnew;
end
f = sum( w(:) .* rho(:) ) + H0;