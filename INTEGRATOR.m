MCI(@(x)abs(x),1000)
T0(@(x)abs(x),10)
function Q=MCI(fun,N)
t=rand(N,1);
t=2*t-1;
f=fun(t);
Q=(2/N)*sum(fun(t));
end
function T=T0(fun,k)
h=2^(2-k);
N=2^(k-1);
T=.5*(fun(-1)+fun(1));
t=fun(linspace(-1+h,1-h,N-1));
T=T+sum(t);
T=h*T;
end