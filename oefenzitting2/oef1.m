%% Error estimation

%% Testfunctions
f1 = @(x) x.^20;
f2 = @(x) exp(x);
f3 = @(x) exp(-x.^2);
f4 = @(x) 1./(1+16*x.^2);

%% Analytische oplossingen
I1 = @(x) x.^21/21;
I2 = @(x) exp(x);
I3 = @(x) (1/2)*sqrt(pi)*erf(x);
I4 = @(x) (1/4)*atan(4*x);

%% Calculate integrals numerically
It1 = ctrap(f1, 31); Ia1 = I1(1)-I1(-1);
It2 = ctrap(f1, 31); Ia2 = I2(1)-I2(-1);
It3 = ctrap(f1, 31); Ia3 = I3(1)-I3(-1);
It4 = ctrap(f1, 31); Ia4 = I4(1)-I4(-1);

fprintf('rel error 1: %f\n', abs(It1 -Ia1)/abs(Ia1));
fprintf('rel error 2: %f\n', abs(It2 -Ia2)/abs(Ia2));
fprintf('rel error 3: %f\n', abs(It3 -Ia3)/abs(Ia3));
fprintf('rel error 4: %f\n', abs(It4 -Ia4)/abs(Ia4));
