function [y] = lattice_reduction_quadratic(x, param)
    % Use the lattice reduction technique in 
    % Bach, Francis. "Submodular functions: from discrete to continuous domains." 
    % Mathematical Programming 175 (2019), Section 4.4
    % to reduce a quadratic 0.5 x^T Q x defined on a lattice with
    % q_ij <= 0, i != j, to a submodular set function.
    u = reshape(x, param.n, param.k-1);
    x_proj = map_row_noninc(u);
    y = lattice_quadratic(sum(x_proj, 2), param) + param.L * param.m * sum(abs(u - x_proj), "all");
end