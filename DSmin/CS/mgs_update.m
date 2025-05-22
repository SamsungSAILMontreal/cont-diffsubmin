function [q_new] = mgs_update(Q, x)
    % Given a matrix Q with orthonormal columns, use modified
    % Gram-Schmidt to orthogonalize x wrt Q. If Q is empty, we just return
    % the normalized x. This code DOES NOT check if Q is orthonormal. 
    if isempty(Q)
        q_new = x / norm(x);
    else
        [~, n] = size(Q);
        v = x;
        for i = 1:n
            v = v - (Q(:, i)' * v) * Q(:, i);
            %atom_new = atom_new - (A_T(:,j)'*atom_new)*A_T(:,j);
        end
        q_new = v / norm(v);
    end
end

