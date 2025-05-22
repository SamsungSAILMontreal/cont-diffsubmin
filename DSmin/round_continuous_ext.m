function [xmin,Fmin] = round_continuous_ext(rho, y0, F_add)
% finds the minimizer of F among the level sets
% (assumes all variables have the same cardinality)


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