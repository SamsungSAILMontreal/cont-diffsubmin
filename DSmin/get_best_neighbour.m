function [x_best, F_best, difference] = get_best_neighbour(F_add, F_rmv, y, x, k, return_self)
% Computes the best neighbouring point of F at x in {0,...,k-1}^n.

% INPUT:
% F_add: Given initial y = F(x), F_add computes y_add = F(x + e_i).
% F_rmv: Given initial y = F(x), F_rmv computes y_add = F(x - e_i).
% x: A point in {0,...,k-1}^n.
% return_self: Use return_self=true if you want to consider x as a
% neighbouring point (default), setting return_self=false will exclude
% x from consideration when choosing a best neighbour. 

% OUTPUT:
% x_best: The best neighbouring point. 
% F_best: F(x_best).
% difference: Array of length 2 with difference(1) giving the corrdinate
% j such that x_best = x +/- e_j. difference(2)=1 means x_best = x + e_j
% while difference(2)=-1 means x_best = x - e_j.       
    if nargin == 5
        return_self = true;
    end
    n = length(x);
    x_best = x;
    F_best = y;
    difference = [0, 0];
    if ~return_self
        F_best = inf;
    end
        for j = 1:n
            if x(j) ~= k-1
                [y_add, x_add] =  F_add(y, x, j);
                if y_add < F_best 
                    F_best = y_add; 
                    x_best = x_add;
                    difference = [j ,1];
                end
            end
    
            if x(j) ~= 0
                [y_rmv, x_rmv] =  F_rmv(y, x, j);
                if y_rmv < F_best
                    F_best = y_rmv; 
                    x_best = x_rmv;
                    difference = [j ,-1];
                end
            end
        end
end