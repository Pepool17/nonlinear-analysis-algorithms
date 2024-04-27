def evaluar_polinomio(a,z):
    z=complex(z)
    n=len(a)-1
    b=[]
    b.append(complex(a[0]))
    for i in range(1,n+1):
        b.insert(i,a[i]+b[i-1]*z)
    resultado=b[n]
    return([resultado,b])


