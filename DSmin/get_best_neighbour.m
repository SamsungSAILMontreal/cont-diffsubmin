function [x_best, F_best, difference] = get_best_neighbour(F_add, F_rmv, y, x, k, return_self)
    % set current sol of F to x_new, for faster computation of marginals
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