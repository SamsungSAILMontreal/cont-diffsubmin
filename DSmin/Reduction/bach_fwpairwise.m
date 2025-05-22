function [rho, gaps] = bach_fwpairwise(F, rho, param, maxiter, eps)
    [n, k] = size(rho);
    k = k+1;
    F0 = F(zeros(n, 1), param);
    w = bach_greedy_algorithm(rho,F,param);

    ws = zeros( n*(k-1) ,maxiter+1 );
    ws(:,1) = reshape(w,n*(k-1),1);
    convex_combinations = zeros(1,maxiter+1);
    convex_combinations(1) = 1;
    iter = 1;
    gaps = inf;
    
    while iter <= maxiter && gaps(end) >= eps

        % compute gradient direction
        rho = zeros(n,k-1);
        for i=1:n
            rho(i,:) = -pav(w(i,:));
        end
        
        % linear oracle
        [wbar,f,Fmin] = bach_greedy_algorithm(rho,F,param);
        ws(:,iter+1) = reshape(wbar,n*(k-1),1);
        % compute away step
        ind = find(convex_combinations(1:iter)>0);
        [a,b] = min( ws(:,ind)'* reshape(rho,n*(k-1),1) );
        b = ind(b);
        away_direction = w - reshape(ws(:,b),n,k-1);
        max_step_away = convex_combinations(b);
        fw_direction = wbar - w;
        max_step_fw = 1;
        direction = fw_direction + away_direction;
        max_step = max_step_away;
        
        
        dual_fw_pair(iter) = .5 * sum( rho(:).^2 ) + sum( rho(:) .* w(:) );
        primal_fw_pair(iter) = f + .5 * sum( rho(:).^2 );
        primal_fw_pair_min(iter) = Fmin;
        dual_fw_pair_min(iter) = sum( min(min( cumsum(w,2) , [],2),0) );
        
        % line search
        aa = sum( rho(:) .* ( direction(:) ) );
        bb = sum( ( direction(:) ).^2 );
        step = min(max_step,max(aa/bb,0));
        convex_combination_direction = zeros(1,maxiter+1);
        convex_combination_direction(iter+1)=1;
        convex_combination_direction(b) = -1;
        convex_combinations = convex_combinations + step * convex_combination_direction;
        w = w + step * direction;

        gaps(iter) = primal_fw_pair_min(iter) - dual_fw_pair_min(iter) - F0;

        iter = iter + 1;
    end
end

