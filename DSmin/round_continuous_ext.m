function [xmin,Fmin] = round_continuous_ext(rho, y0, F_add)
% Given an n*(k-1) row-nonincreasing matrix rho with elements in [0,1],
% round rho to a point xmin in {0,...,k-1}^n such that F(xmin) is less than
% or equal to the continuous extension of F at rho.

% INPUT:
% rho: An n*(k-1) row-nonincreasing matrix with elements in [0,1].
% y0: y0 = F(0).
% F_add: Given initial y = F(x), F_add computes y_add = F(x + e_i).

% OUTPUT:
% xmin: The rounded point in {0,...,k-1}^n such that F(xmin) is less than
% or equal to the continuous extension of F at rho.
% Fmin: F(xmin).


[n, k] = size(rho);
k = k + 1;

% first order all rhos (does preserve the ordering within rows if equal
% values)
[~, s] = sort( reshape(rho',n*(k-1),1), 'descend');
[js,is] = ind2sub([k-1 n],s);

% now go through all elements
xold = zeros(n,1);
Fold = y0;

xmin = xold;
Fmin = Fold;

for i=1:n*(k-1)
    xnew = xold; 
    xnew(is(i))=js(i);
    Fnew = F_add(Fold, xold, is(i));
    if (Fnew<Fmin) 
        xmin = xnew; 
        Fmin = Fnew; 
    end
    xold = xnew;
    Fold = Fnew;
end