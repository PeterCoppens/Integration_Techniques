function q = apply_rule(w, f)
    % apply quadrature rule or null rule to f where f is a function handle
    % and w is the row vector representation of the rule.
    k = length(w)-1;
    fx = f(-1 + (2/k)*(0:k));
    q = w*fx';
end

