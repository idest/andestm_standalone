import numpy as np

f1 = 0.5
f2 = 0.5

n1 = 3 
A1 = 7.00e-14
H1 = 520
logA1 = 4.5

n2 = 3.8
A2 = 4.5e-38
H2 = 8.9
logA2 = -6.9

def poly_flow(f1,f2,n1,A1,H1,n2,A2,H2):
    n = 10**(f1*np.log10(n1)+f2*np.log10(n2))
    H = (H2*(n-n1) - H1*(n-n2))/(n2-n1)
    A = 10**((np.log10(A2)*(n-n1)-np.log10(A1)*(n-n2))/(n2-n1))
    return [n, H, A]

def poly_flow_n(f1,f2,n1,A1,H1,n2,A2,H2):
    n = np.exp(f1*np.log(n1)+f2*np.log(n2))
    H = (H2*(n-n1) - H1*(n-n2))/(n2-n1)
    A = np.exp((np.log(A2)*(n-n1)-np.log(A1)*(n-n2))/(n2-n1))
    return [n, H, A]

def poly_flow_log(f1,f2,n1,logA1,H1,n2,logA2,H2):
    n = 10**(f1*np.log10(n1)+f2*np.log10(n2))
    H = (H2*(n-n1) - H1*(n-n2))/(n2-n1)
    logA = ((logA2)*(n-n1)-(logA1)*(n-n2))/(n2-n1)
    return [n, H, logA]

print(poly_flow(f1,f2,n1,A1,H1,n2,A2,H2))
