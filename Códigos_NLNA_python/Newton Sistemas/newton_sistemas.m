function [x,nit]=newton_sistemas(F,J,x0,toll,nmax,p)

n=length(F);
nit=0;
Fxn=zeros(n,1);
x=x0;
err=toll+1;

for i=1:n
    for j=1:n
        Jxn(i,j)=eval(J(i-1)*n+j,:);
    end
end

[L,U,P]=lu(Jxn);
step=0;

while err>toll
    if step==p
        step=0;
        for i=1:n
            Fxn(i)=eval(F(i,:));
            for j=1:n
                Jxn(i,j)=eval(J((i-1)*n+j,:));
            end
        end
        [L,U,P]=lu(Jxn);
    else
        for i=1:n
            Fxn(i)=eval(F(i,:));
        end
    end
    nit=nit+1;
    step=step+1;
    Fxn=-P*Fxn;
    y=forward_col(L,Fxn);
    deltax=backward_col(U,y);
    x=x+deltax;
    err=norm(deltax);
    if nit>nmax
        disp('Se alcanzo el numero maximo de it sin conv');
        break
    end
end