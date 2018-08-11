from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.misc.misc_c import prod
from sage.matrix.constructor import matrix
def colloc(c,P):
    """
    Given a list C of collocation points, and a ring R (AA), build 
    the A and B part of the Butcher array of an associated  Runge-Kutta
    method.

    AUTHOR::

    Thierry Dumont (2016).

    EXAMPLES::

    sage: R = PolynomialRing(AA,"x")
    sage: n = 4
    sage: c = [(s[0]+1)/2 for s in R(legendre_P(n,x)).roots()]
    sage: A,B = colloc(c,R)
    """
    Pb= P.base()
    x = P.gen()
    n = len(c)
    #
    for s in c:
        if s<0 or s>1:
            raise IndexError("colloc: collocation point ",s," is out of [0,1]")
    pols=[]
    for  i in range(0,n):
        ploc=P(1)
        for p in range(0,n):
            if i!=p:
                ploc*=P((x-c[p])/(Pb(c[i])-Pb(c[p])))
        pols.append(ploc)

    prims = [p.integral(x) for p in pols]
    prims0 = [p(x = 0) for p in prims]
    A = matrix(Pb,[ [prims[j](x = c[i])-prims0[j] for j in range(0,n)] \
                  for i in range(0,n)])
    
    B = [prims[j](x = 1) - prims0[j] for j in range(0,n)]
    # exactify to improve lisibility, if possible!
    for i in range(0,n):
        for j in range(0,n):
            A[i,j].exactify()
        B[i].exactify()
    return A,B

