
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
n=15;
y=linspace(0,1,n+1);
x=linspace(-1,1,2*n+1);
r=reductions(f1,x)
figure(1)
subplot(2,2,1)
semilogy(1:length(errors(f1,x)),abs(errors(f1,x)),'--o')
title('ej for x^{20}')
subplot(2,2,2)
semilogy(1:length(errors(f2,x)),abs(errors(f2,x)),'--o')
title('ej for exp(x)')
subplot(2,2,3)
semilogy(1:length(errors(f3,x)),abs(errors(f3,x)),'--o')
title('ej for exp(-x^2)')
subplot(2,2,4)
semilogy(1:length(errors(f4,x)),abs(errors(f4,x)),'--o')
title('ej for 1/(1+16x^2)')
t=0;
figure(2)
subplot(2,2,1)
plot(1:length(reductions(f1,x)),reductions(f1,x),'--o')
title('r_j for x^{20}')
subplot(2,2,2)
plot(1:length(reductions(f2,x)),reductions(f2,x),'--o')
title('r_j for exp(x)')
subplot(2,2,3)
plot(1:length(reductions(f3,x)),reductions(f3,x),'--o')
title('r_j for exp(-x^2)')
subplot(2,2,4)
plot(1:length(reductions(f4,x)),reductions(f4,x),'--o')
title('r_j for 1/(1+16x^2)')
%adapted reduction factors
figure(3)
subplot(2,2,1)
plot(1:length(reductions2(f1,x)),reductions2(f1,x),'--o')
title('R_j for x^20')
subplot(2,2,2)
plot(1:length(reductions2(f2,x)),reductions2(f2,x),'--o')
title('R_j for exp(x)')
subplot(2,2,3)
plot(1:length(reductions2(f3,x)),reductions2(f3,x),'--o')
title('R_j for exp(-x^2)')
subplot(2,2,4)
plot(1:length(reductions2(f4,x)),reductions2(f4,x),'--o')
title('R_j for 1/(1+16x^2)')
U=nullrules(y);
y2=fliplr(y)
y=[y2(1:end-1),y];
U2=fliplr(U);
U=[U2(1:end,1:end-1),U];
I=U*(x)'
%% Calculate Null rules
% Vandermonde matrix
function [ U ] = nullrules( x )
    n = length(x); % n+1 in fact
    U = zeros(n);
    vandermonde = flipud(vander(x)');
    vandermonde=vandermonde.^2;
    vandermonde(:,1)=0;
    Vandermonde(1,1)=.5;
    for m = 1:n
        vandermonde = vandermonde(1:end-1,:); 
        nullSpace = null(vandermonde);
        nullVector = nullSpace(:,1);
        for i = 1:m-1
            factor = nullVector'*U(:,i);
            nullVector = nullVector - factor*U(:,i);
        end
        U(:,m) = nullVector;
        U(:,m) = U(:,m)/norm(U(:,m));
    end
    U=flipud(U');
    U=U(2:end-1,1:end);
end
function I=CheckNullRule(U,x)
I=zeros(size(U,1),1);
for i=1:size(I)
    I(i)=U(i,:)*(x.^(i-1))';
end
end
function e=errors(f,x)
U=nullrules(x((end+1)/2:end));
U2=fliplr(U);
U=[U2(1:end,1:end-1),U];
e=U*feval(f,x)';
ind=find(abs(e)>10^(-14));
e=e(ind);
end
function E=errorsnew(f,x)
e=errors(f,x);
n=(length(e)-mod(length(e),2) )/2
for i=1:n
    E(i)=sqrt(e(2*i-1)^2+e(2*i)^2);
end
end

function r=reductions(f,x)
e=errors(f,x);
r=e(2:end)./e(1:end-1);
end
function R=reductions2(f,x)
E=errorsnew(f,x);
R=E(2:end)./E(1:end-1);
end

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
