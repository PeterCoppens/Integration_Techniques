%% Gauss-Legendre, Clenshaw-Curtis and Romberg

%% Testfunctions
f1 = @(x) x.^20;
f2 = @(x) exp(x);
f3 = @(x) exp(-x.^2);
f4 = @(x) 1./(1+16*x.^2);
f5 = @(x) exp(-x.^(-2));
f6 = @(x) abs(x).^3;

%% Plots
x = linspace(-1, 1, 100);
figure(); hold on;
plot(x, f1(x));
plot(x, f2(x));
plot(x, f3(x));
plot(x, f4(x));
plot(x, f5(x));
plot(x, f6(x));

%% Analytische oplossingen
I1 = @(x) x.^21/21;
I2 = @(x) e.^x;
I3 = @(x) (1/2)*sqrt(pi)*erf(x);
I4 = @(x) (1/4)*atan(4*x);
I5 = @(x) sqrt(pi)*erf(1./x) + exp(-1./x.^2).*x;
I6 = @(x) (1/4).*x.^4*sgn(x);

