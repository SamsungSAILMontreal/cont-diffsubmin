function [A] = map_row_noninc(A)
    % Map an arbitrary m by n binary matrix A 
    % onto the set {0,1}↓^{m x n} according to the rule in 
    % Bach, Francis. "Submodular functions: from discrete to continuous domains." 
    % Mathematical Programming 175 (2019) 
    % Section 4.4
    for i = 1:size(A, 1)
        j = find(A(i, :), 1, 'last');
        if ~isempty(j)
            A(i, 1:j) = 1;
        end
    end
end