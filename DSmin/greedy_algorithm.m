function [w, f, Fmin, xmin] = greedy_algorithm(rho, y0, F_add, ties)
% greedy algorithm (assumes all variables have the same cardinality)

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