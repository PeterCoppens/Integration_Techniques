function q = ctrap(k, f)    
    h = 2/k;
    q = (h/2)*[1, 2*ones(1, k-2), 1];
    
    if exist('f','var')
        x = -1 + h*(0:k);
        xf = linspace(-1, 1, 1000);
        plot(xf, abs(f(xf)-interp1(x,f(x),xf)), 'b-');
    end
end

