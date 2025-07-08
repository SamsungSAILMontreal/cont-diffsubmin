function [y_new, x_new] = submod_major_marginal(Q, c, lambda, W, y, x, grid, i, direction)
% Computes the marginal gain/loss of
% f(x) = 0.5*grid(x + 1)*Q*grid(x + 1)' + c^T*grid(x + 1)' + lambda * ||x||_0 + mod_F(x)
% with mod_F a modular function characterized by W. 
% Given input/value pair y = f(x), compute y_new = f(x +/- e_i), with
% e_i the ith basis vector.

% INPUT:
% W: An n*(k-1) matrix.
% y: The current value y = f(x).
% x: A point in {0,...,k-1}^n.
% grid: row vector of length (k-1) containing the (ordered) grid points.
% direction: direction="add" computes y_new = f(x + e_i), 
%            direction="rmv" computes y_new = f(x - e_i).

% OUTPUT:
% x_new: x_new = x + e_i if direction="add", else x_new = x - e_i if
% direction="rmv".
% y_new: F(x_new).
    if strcmp(direction, "add")
        if x(i) < length(grid)
            x_new = x;
            x_new(i) = x_new(i) + 1;
            w = W(i, x_new(i));
            d = grid(x(i) + 2) - grid(x(i) + 1);
            y_new = y + 0.5*d*grid(x + 1)*(Q(:,i)+Q(i,:)') + 0.5*d^2 * Q(i,i) + d*c(i) + w;
            if grid(x(i) + 2) == 0
                y_new = y_new - lambda;
            end
            if grid(x(i) + 1) == 0 && grid(x(i) + 2) ~= 0
                y_new = y_new + lambda;
            end
        else
            x_new = x;
            y_new = y;
        end
    elseif strcmp(direction, "rmv")
        if x(i) > 0
            w = W(i, x(i));
            x_new = x;
            x_new(i) = x_new(i) - 1;
            d = grid(x(i) + 1) - grid(x(i));
            y_new = y - 0.5*d*grid(x + 1)*(Q(:,i)+Q(i,:)') + 0.5*d^2 * Q(i,i) - d*c(i) - w;
            if grid(x(i)) == 0
                y_new = y_new - lambda;
            end
            if grid(x(i) + 1) == 0 && grid(x(i)) ~= 0
                y_new = y_new + lambda;
            end
        else
            x_new = x;
            y_new = y;
        end
    else
        error("Invalid Direction")
    end
end

