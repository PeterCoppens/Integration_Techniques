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
I2 = @(x) exp(x);
I3 = @(x) (1/2)*sqrt(pi)*erf(x);
I4 = @(x) (1/4)*atan(4*x);
I5 = @(x) sqrt(pi)*erf(1./x) + exp(-1./(x.^2)).*x;
I6 = @(x) (1/4).*x.^4*sign(x);

%% Calculate integrals numerically
n = 1:100;
errg = zeros(6, length(n));
errcc = zeros(6, length(n));
for i = n
    Iex = I1(1)-I1(-1);
    Ig = gauss(f1, i);
    Icc = clenshaw_curtis(f1, i);
    errg(1, i) = abs(Ig-Iex)/abs(Iex);
    errcc(1, i) = abs(Icc-Iex)/abs(Iex);
    
    Iex = I2(1)-I2(-1);
    Ig = gauss(f2, i);
    Icc = clenshaw_curtis(f2, i);
    errg(2, i) = abs(Ig-Iex)/abs(Iex);
    errcc(2, i) = abs(Icc-Iex)/abs(Iex);
    
    Iex = I3(1)-I3(-1);
    Ig = gauss(f3, i);
    Icc = clenshaw_curtis(f3, i);
    errg(3, i) = abs(Ig-Iex)/abs(Iex);
    errcc(3, i) = abs(Icc-Iex)/abs(Iex);
    
    Iex = I4(1)-I4(-1);
    Ig = gauss(f4, i);
    Icc = clenshaw_curtis(f4, i);
    errg(4, i) = abs(Ig-Iex)/abs(Iex);
    errcc(4, i) = abs(Icc-Iex)/abs(Iex);
    
    Iex = I5(1)-I5(-1);
    Ig = gauss(f5, i);
    Icc = clenshaw_curtis(f5, i);
    errg(5, i) = abs(Ig-Iex)/abs(Iex);
    errcc(5, i) = abs(Icc-Iex)/abs(Iex);
    
    Iex = I5(1)-I5(-1);
    Ig = gauss(f5, i);
    Icc = clenshaw_curtis(f5, i);
    errg(5, i) = abs(Ig-Iex)/abs(Iex);
    errcc(5, i) = abs(Icc-Iex)/abs(Iex);
    
    Iex = I6(1)-I6(-1);
    Ig = gauss(f6, i);
    Icc = clenshaw_curtis(f6, i);
    errg(6, i) = abs(Ig-Iex)/abs(Iex);
    errcc(6, i) = abs(Icc-Iex)/abs(Iex);
end

%% Plot errors
figure();
for i=1:6
    subplot(3,2,i); semilogy(errg(i, :)); hold on;
    subplot(3,2,i); semilogy(errcc(i, :));
end

%% Plot errors
figure();
fprintf('amount of points required to get 7 correct digits:\n');
fprintf('==================================================\n');
for i=1:6
    subplot(3,2,i); plot(-log10(errg(i, :))); hold on;
    subplot(3,2,i); plot(-log10(errcc(i, :)));
    cost_g = find(-log10(errg(i, :))>=7, 1, 'first');
    cost_cc = find(-log10(errg(i, :))>=7, 1, 'first');
    fprintf('f%d: cost -- gauss: %d -- clenshaw-curtis: %d\n', i, n(cost_g), n(cost_cc));
end
