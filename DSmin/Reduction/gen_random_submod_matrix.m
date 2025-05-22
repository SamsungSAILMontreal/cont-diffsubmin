function [A] = gen_random_submod_matrix(n)
    % Generate a random symmetric n by n matrix A with the properties that:
    % 1. a_ij ~ U(-1/n, 0) for i != j.
    % 2. a_ii ~ U(0, 1).
    M = -rand(n) / n;
    upper_M = triu(M, 1);
    A = upper_M + upper_M' + diag(rand(1,n),0);
end