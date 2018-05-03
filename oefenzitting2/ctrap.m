function Q = ctrap(f, k)
    a = -1; b = 1;
    xj = [a, b];
    xi = linspace(a, b, k+1);
    xij = @(i, j) xi(i-1) + ((b-a)/(2*k))*(1+xj(j));
    w = [(b-a)/2, (b-a)/2];

    Q = 0;
    for i=2:k+1
        for j=1:2
            Q = Q + w(j)*f(xij(i, j));
        end
    end
    Q = ((b-a)/(2*k))*Q;
end

