MCI(@(x)abs(x),1000)
T0(@(x)abs(x),5)
T_rec(@(x)abs(x),1,4)
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
function T=T_rec(fun,j,m)
if (j==0)
    T=T0(fun,m);
else
    T=(4^j*T_rec(fun,j-1,m+1)-T_rec(fun,j-1,m))/((4^j)-1);
end
end
