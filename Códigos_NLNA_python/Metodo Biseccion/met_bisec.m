function [vptmed,verr,x,fx,nit,err]=met_bisec(a,b,tolle,tollc,nmax,fun) 

f=inline(fun);
fa=f(a);
fb=f(b);

c=(a+b)/2;
fc=f(c);

vptmed=[];  
verr=[];
err=tolle+1; 
nit=0;

if (sign(fa)==sign(fb)) 
    ezplot(f,[a,b]);
    error('no hay cero de f en [a,b]'); 
elseif abs(fc)==0
    ezplot(f,[a,b]);
    fprintf('(a+b)/2 es una raiz o esta cerca de una raiz')
    x=c; fx=fc; verr=abs(fc); err=abs(fc); nit=0; vptmed=0;
    return
elseif abs(fa)==0
    ezplot(f,[a,b]);
    fprintf('a es una raiz o esta cerca de una raiz')
    x=a; fx=fa; verr=abs(fa); err=abs(fa); nit=0; vptmed=0;
    return
elseif abs(fb)==0
    ezplot(f,[a,b]);
    fprintf('b es una raiz o esta cerca de una raiz')
    x=b; fx=fb; verr=abs(fb); err=abs(fb); nit=0; vptmed=0;
    return
end
 

while (nit < nmax && err > tolle && abs(fc)>tollc)
     nit=nit+1; 
     fa=f(a);
     fb=f(b);
     e=(b-a)/2;
     c=a+e;
     fc=f(c);
     vptmed=[vptmed;c]; 

     if (sign(fa)~=sign(fc)) 
         b=c; 
     elseif (sign(fc)~=sign(fb)) 
         a=c;
     end
     
     err=abs(b-a);
     verr=[verr;err];   
end
x=c;
fx=fc;

figure(1)
ezplot(f,[c-1,c+1])

figure(2)
semilogy(verr)

figure(3)
plot(vptmed)