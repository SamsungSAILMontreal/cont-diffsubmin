function [rho, gaps] = fwpairwise(y0, F_add, rho, maxiter, eps)
%A modified version of the Pairwise Frank-Wolfe algorithm for submodular 
% minimization from "Submodular Functions: from Discrete to Continous
% Domains". Use for minimizing submodular functions
% min_{x in {0,...,k-1}^n} F(x).

% INPUT:
% y0: y0 = F(0).
% F_add: A function capable of computing the marginal gain of F. That is,
%        given initial y = F(x), F_add computes y_add = F(x + e_i).
% rho: An n*(k-1) row-nonincreasing matrix with elements in [0,1] which 
%      serves as a starting point for Pairwise FW.
% maxiter: max number of Pairwise FW iterations.

% OUTPUT:
    [n, k] = size(rho);
    k = k+1;
    w = greedy_algorithm(rho, y0, F_add);
    F0 = y0;

    ws = zeros( n*(k-1) ,maxiter+1 );
    ws(:,1) = reshape(w,n*(k-1),1);
    convex_combinations = zeros(1,maxiter+1);
    convex_combinations(1) = 1;
    gaps = zeros(1, maxiter);
    
    for iter=1:maxiter

        % compute gradient direction
        rho = zeros(n,k-1);
        for i=1:n
            rho(i,:) = -pav(w(i,:));
        end
        
        % linear oracle
        [wbar, f, Fmin] = greedy_algorithm(rho, y0, F_add);
        ws(:,iter+1) = reshape(wbar,n*(k-1),1);
        % compute away step
        ind = find(convex_combinations(1:iter)>0);
        [a,b] = min( ws(:,ind)'* reshape(rho,n*(k-1),1) );
        b = ind(b);
        away_direction = w - reshape(ws(:,b),n,k-1);
        max_step_away = convex_combinations(b);
        fw_direction = wbar - w;
        direction = fw_direction + away_direction;
        max_step = max_step_away;
        
        primal_fw_pair_min = Fmin;
        dual_fw_pair_min = sum( min(min( cumsum(w,2) , [],2),0) );
        
        % line search
        aa = sum( rho(:) .* ( direction(:) ) );
        bb = sum( ( direction(:) ).^2 );
        step = min(max_step,max(aa/bb,0));
        convex_combination_direction = zeros(1,maxiter+1);
        convex_combination_direction(iter+1)=1;
        convex_combination_direction(b) = -1;
        convex_combinations = convex_combinations + step * convex_combination_direction;
        w = w + step * direction;
        gap = primal_fw_pair_min - dual_fw_pair_min - F0;
        gaps(iter) = gap;

        if gap < eps
            break
        end

    end
    gaps = gaps(1:iter);
end