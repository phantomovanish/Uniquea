"""Numerical check: correlations that are a positive-definite function of the angle between
measurement directions (an isotropic plenum) never exceed |S| = 2*sqrt(2); local (deterministic) models
never exceed 2. Usage: python3 check.py"""
import numpy as np
from numpy.polynomial import legendre as L
from scipy.optimize import minimize
rng=np.random.default_rng(1)
def unit(t,p): return np.array([np.sin(t)*np.cos(p),np.sin(t)*np.sin(p),np.cos(t)])
def S_of(x,K):
    a,a2,b,b2=[unit(x[2*i],x[2*i+1]) for i in range(4)]
    E=lambda u,v: K(np.clip(u@v,-1,1))
    return E(a,b)+E(a,b2)+E(a2,b)-E(a2,b2)
def maxS(K,starts=40):
    best=0
    for _ in range(starts):
        r=minimize(lambda x:-abs(S_of(x,K)),rng.uniform(0,np.pi,8),method="Nelder-Mead",options={"maxiter":4000,"xatol":1e-9,"fatol":1e-12})
        best=max(best,-r.fun)
    return best
print("Kernel E = cos(angle) (spin-1/2 singlet up to sign):  max|S| = %.6f   2*sqrt(2) = %.6f"%(maxS(lambda c:c),2*np.sqrt(2)))
worst=0
for k in range(25):
    c=rng.dirichlet(np.ones(6))           # non-negative Legendre weights, sum 1 -> positive-definite kernel, K(0)=1
    K=lambda x,c=c: L.legval(x,np.r_[0,c] if False else c)
    m=maxS(K,15); worst=max(worst,m)
print("25 random positive-definite isotropic kernels (Legendre P0..P5, weights >= 0, sum 1): largest max|S| = %.6f"%worst)
# local deterministic: outcomes +-1 fixed per setting -> enumerate 16 strategies
best=max(abs(A*B+A*B2+A2*B-A2*B2) for A in(1,-1) for A2 in(1,-1) for B in(1,-1) for B2 in(1,-1))
print("Local deterministic strategies: max|S| = %d"%best)
