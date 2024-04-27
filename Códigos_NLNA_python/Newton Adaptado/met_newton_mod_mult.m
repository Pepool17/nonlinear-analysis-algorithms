function [xvect,fx,x,nit,err,m]=met_newton_mod_mult(fun,dfun,x0,nmax,toll)

f=inline(fun);
df=inline(dfun);
xvect=x0; 
nit=0; 
r=(1); 
m=(1); 
xdif=[];

err=toll+1; 

while (nit < nmax && err > toll)
    nit=nit+1; 
    x=xvect(nit); 
    fx(nit)=f(x);%eval(fun); 
    dfx=df(x);%eval(dfun);
    if (dfx == 0)
        disp('Stop due to vanishing derivative'); 
        return; 
    end
    x=x-m(nit)*fx(nit)/dfx; 
    xvect=[xvect;x]; 
    fx=[fx;eval(fun)];
    rd=err; 
    err=abs(xvect(nit+1)-xvect(nit)); 
    xdif=[xdif;err];
    ra=err/rd; 
    r=[r;ra]; 
    diff=abs(r(nit+1)-r(nit));
    if (diff < 1.e-3 && r(nit+1) > 1.e-2)
        m(nit+1)=max(m(nit),1/abs(1-r(nit+1))); 
    else
        m(nit+1)=m(nit); 
    end
end

if nit>nmax
    ezplot(f,[x0-1,x0+1])
    fprintf('el metodo no converge a la tolerancia deseada y se alcanzo el numero maximo de it\n');
end

figure(1)
ezplot(f,[x-1,x+1])

figure(2)
semilogy(xdif)

figure(3)
plot(xvect)
