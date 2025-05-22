function [y] = lattice_quadratic(x, param)
    % Evaluate the quadratic 0.5 * v' * param.Q * v with v_i's in [lb, ub].
    % The discretization of [lb, ub] is handled by param.m, with
    % m = (ub-lb)/(k-1).
    v = param.lb + param.m * x;
    y = 0.5 * v' * param.Q * v;
end