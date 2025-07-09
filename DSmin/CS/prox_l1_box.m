function out = prox_l1_box(x,alpha)
    %computes the proximal operator of the function 
    %alpha*||x||_1 + delta_{[-1,1]}(x)
    %where delta_{[-1,1]^n}(x) is the convex indicator function of 
    %[-1,1]^n.
    out = max(min(abs(x) - alpha, 1), 0).* sign(x);
end

