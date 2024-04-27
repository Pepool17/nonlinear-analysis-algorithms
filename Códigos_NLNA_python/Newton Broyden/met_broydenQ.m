function [x,it]=met_broydenQ(x,Q,nmax,toll,f)
n=length(f);
it=0;
err=1;

fk=zeros(n,1);
fk1=fk;
for i=1:n
    fk(i)=eval(f(i,:));
end

while it<nmax && err>toll
    s=-Q\fk;
    x=x+s;
    err=norm(s,inf);
    if err>toll
        for i=1:n
            fk1(i)=eval(f(i,:));
        end
        Q=Q+(1/(s'*s))*fk1*s';
    end
    it=it+1;
    fk=fk1;
end