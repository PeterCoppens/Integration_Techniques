%% Testfunctions
f1 = @(x) x.^20;
f2 = @(x) exp(x);
f3 = @(x) exp(-x.^2);
f4 = @(x) 1./(1+16*x.^2);
f = {f1, f2, f3, f4};

%% Error estimate
n = [-0.5, 0.5];
q = [0.5, 0.5];
xh = linspace(-1, 1, 31);
h = xh(2)-xh(1);
Et = zeros(4, 1);
for j = 1:4
    for i = 1:length(xh)-1
        fh = @(x) f{j}(0.5*(x*(xh(i)+xh(i+1)) + xh(i+1)-xh(i)));
        Ih = @(x) I{j}(0.5*(x*(xh(i)+xh(i+1)) + xh(i+1)-xh(i)));
        e1 = apply_rule(n, fh)*0.5*h;
        ft = apply_rule(q, fh);
        et = apply_rule(q, @(x) abs(fh(x) - ft))*0.5*h;
        r1 = abs(e1)/et;
        if r1 > 1/200
            Et(j) = Et(j) + et;
        else
            Et(j) = Et(j) + (200^(1.5))*(r1^(1/2))*abs(e1);
        end
    end
end

disp(Et);

function q = apply_rule(w, f)
    % apply quadrature rule or null rule to f where f is a function handle
    % and w is the row vector representation of the rule.
    k = length(w)-1;
    fx = f(-1 + (2/k)*(0:k));
    q = w*fx';
end

