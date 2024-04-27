from Evaluar_polinomio import evaluar_polinomio
import numpy as np

def newton_horner(a,x0,tol,max_it):
    x0=complex(x0)
    n=len(a)-1
    raices=[]
    for _ in range(n):
        it=0
        err=1

        while err>tol and it<max_it:
            fx = evaluar_polinomio(a,x0)
            if fx[0]==0:
                    err=0
                    it += 1
            dfx = evaluar_polinomio(fx[1][:-1],x0)
            xn = x0 - fx[0]/dfx[0]
            err = abs(xn-x0)
            x0 = xn   
            it += 1
        raices.append(xn)
        fx = evaluar_polinomio(a,x0)
        a = fx[1][:-1]
        if xn==0:
            x0=1+1j
    return raices

if __name__=="__main__":
    a=[1,0,1]
    tol=1e-20
    x0=2+1j
    max_it=1000
    raices=newton_horner(a,x0,tol,max_it)
    
    print("raices")
    for raiz in raices:
        print(raiz)

