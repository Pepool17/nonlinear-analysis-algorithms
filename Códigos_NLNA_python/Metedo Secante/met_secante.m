function [verr,fx,x,nit,err]=met_secante(a,b,x1,x0,nmax,toll,fun) 

f=inline(fun);
fa=f(a);
fb=f(b); 
m=(b-a)/(fb-fa);

xkm=x0; 
xk=x1;

verr=[]; 
err=toll+1; 
nit=0; 

while (nit < nmax && err > toll)
    nit=nit+1; 
    fxk=f(xk);
    fxkm=f(xkm); 
    x=xk-((xk-xkm)/(fxk-fxkm))*fxk;
    xkm=xk;
    xk=x;
    fxk=f(xk);
    fxkm=f(xkm);
    err=abs(xk-x);
    verr=[verr;err]; 
    x=xn; 
end

figure(1)
ezplot(f,[x-1,x+1])

figure(2)
semilogy(verr)

figure(3)
plot(xvect)