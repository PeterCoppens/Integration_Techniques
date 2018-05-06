%% Error estimation

%% Testfunctions
f1 = @(x) x.^20;
f2 = @(x) exp(x);
f3 = @(x) exp(-x.^2);
f4 = @(x) 1./(1+16*x.^2);
f = {f1, f2, f3, f4};

%% Analytische oplossingen
I1 = @(x) x.^21/21;
I2 = @(x) exp(x);
I3 = @(x) (1/2)*sqrt(pi)*erf(x);
I4 = @(x) (1/4)*atan(4*x);
I = {I1, I2, I3, I4};

%% Calculate quadrature rule and null rules
q = ctrap(30);
U = null_moment(q, floor(length(q)/2));

%% Plot behaviour
figure;
for i=1:4
    subplot(2,2,i); hold on;
    title(sprintf('$f_%d(x)$', i), 'interpreter', 'latex');
    ctrap(30, f{i});
end

%% Calculate integrals
It = zeros(1, 4);
Ia = zeros(1, 4);
for i=1:4
    It(i) = apply_rule(q, f{i}); Ia(i) = I{i}(1)-I{i}(-1);
end

%% Print relative errors
Ea = zeros(4, 1);
fprintf('relative errors: \n========================== \n');
for i=1:4
    Ea(i) = abs(It(i)-Ia(i));
    fprintf('rel error %d: %f\n', i, abs(It(i) -Ia(i))/abs(Ia(i)));
end
fprintf('\nabs errors: \n========================== \n');
for i=1:4
    fprintf('rel error %d: %f\n', i, Ea(i));
end

%% Validate null rules
r = zeros(2, size(U, 1));   % residue
tol = 1e-13;
r(1, :) = 10*tol;

% initialize vandermonde matrix
k = size(U, 1);
x = -1 + (2/k)*(0:k);
V = flipud(vander(x)');

for i=1:size(U, 1)
    d = size(V, 1)+1;
    while r(1,i) > tol && d >= 1
        d = d-1;
        r(1,i) = norm(V(1:d, :)*(U(i, :)'));
    end
    r(2,i) = d;
end

%% Plot the error estimates
E = zeros(4, size(U, 1));
for i=1:4
    E(i, :) = abs(apply_rule(U, f{i}))';
end

figure;
for i = 1:4
    semilogy(r(2, :)', E(i, :)'); hold on;
end
xlabel('$d$', 'interpreter', 'latex');
ylabel('$e_j$', 'interpreter', 'latex');
legend({'$f_1$', '$f_2$', '$f_3$', '$f_4$'}, 'interpreter', 'latex');

%% Reduction factors
R = zeros(size(E)-[0, 1]);
for i = 1:size(E, 2)-1
    check = E(:, i+1)>10^(-10);
    R(:, i) = check.*(E(:,i)./E(:,i+1)) + ~check.*(10^(-10));
end
Rm = max(R, [], 2);

%% Remove phase effect
Ep = zeros(4, floor(size(E, 2)/2));
for i=1:size(Ep, 2)
    Ep(:, i) = sqrt(E(:, 2*i-1).^2+E(:, 2*i).^2);
end

figure();
semilogy(2*(1:size(Ep, 2))-1, flipud(Ep'));
xlabel('$d$', 'interpreter', 'latex');
ylabel('$E_j$', 'interpreter', 'latex');
legend({'$f_1$', '$f_2$', '$f_3$', '$f_4$'}, 'interpreter', 'latex');

%% Reduction factors
Rp = zeros(size(Ep)-[0, 1]);
for i = 1:size(Ep, 2)-1
    check = Ep(:, i+1)>10^(-10);
    Rp(:, i) = check.*(Ep(:,i)./Ep(:,i+1)) + ~check.*(10^(-10));
end
Rpm = max(Rp, [], 2);

%% Plot the reduction factors
figure();
for i=1:4 
    subplot(211); hold on;
    idx = fliplr(1:size(R, 2));
    plot(idx(R(i, :)~=0), R(i, R(i, :)~=0)', '.-');
    subplot(212); hold on;
    idx = fliplr(1:size(Rp, 2));
    plot(2*idx(Rp(i, :)~=0)-1, Rp(i, Rp(i, :)~=0)', '.-');
end
subplot(211);
xlabel('$d$', 'interpreter', 'latex');
ylabel('$r_j$', 'interpreter', 'latex');
legend({'$f_1$', '$f_2$', '$f_3$', '$f_4$'}, 'interpreter', 'latex');
subplot(212);
xlabel('$d$', 'interpreter', 'latex');
ylabel('$R_j$', 'interpreter', 'latex');
legend({'$f_1$', '$f_2$', '$f_3$', '$f_4$'}, 'interpreter', 'latex');

%% Error estimate
n = [-0.5, 0.5];
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