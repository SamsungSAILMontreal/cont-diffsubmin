function out = prox_l1_box(x,alpha)
    out = max(min(abs(x) - alpha, 1), 0).* sign(x);
end

