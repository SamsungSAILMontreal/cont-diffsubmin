function [y_new, x_new] = F_marginal(Q, c, lambda, y, x, grid, i, direction)
    if strcmp(direction, "add")
        if x(i) < length(grid) - 1
            x_new = x;
            x_new(i) = x_new(i) + 1;
            d = grid(x(i) + 2) - grid(x(i) + 1);
            y_new = y + 0.5*d*grid(x + 1)*(Q(:,i)+Q(i,:)') + 0.5*d^2 * Q(i,i) + d*c(i);
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
            x_new = x;
            x_new(i) = x_new(i) - 1;
            d = grid(x(i) + 1) - grid(x(i));
            y_new = y - 0.5*d*grid(x + 1)*(Q(:,i)+Q(i,:)') + 0.5*d^2 * Q(i,i) - d*c(i);
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

