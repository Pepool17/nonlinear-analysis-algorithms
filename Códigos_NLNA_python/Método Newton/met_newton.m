function [xvect,verr,fx,x,nit,err]=met_newton(fun,dfun,x0,M,tol)

f=inline(fun);
df=inline(dfun);

xvect=x0; 
x=x0; 

fx=f(x);
dfx=df(x);

verr=[]; 
err=tol+1;
nit=0;

while (err>tol && nit<M)
    nit=nit+1;
    diff=fx/dfx;
    x=x-diff;
    xvect=[xvect;x];
    diff=abs(diff);
    fx=f(x);
    dfx=df(x);
    err=diff;
    verr=[verr;err];
end

if nit>M
    ezplot(f,[x0-1,x0+1])
    fprintf('el metodo no converge a la tolerancia deseada y se alcanzo el numero maximo de it\n');
end

figure(1)
ezplot(f,[x-1,x+1])

figure(2)
semilogy(verr)

figure(3)
plot(xvect)
