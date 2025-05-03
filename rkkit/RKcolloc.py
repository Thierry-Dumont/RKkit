from sage.all import *
from .RKRungeKutta import RungeKutta
from .RKExceptions import *
def colloc(c,P,title):
    """
    Given a list C of collocation points in [0,1], and an univariate polynomial
    ring P (over an exact field -actually over AA-),  build the A and B part
    of the Butcher array of an associated  Runge-Kutta method and return a
    Runge-Kutta method class.

    The "title" parameter is the name given to the generated Runge-Kurtta
    method.

    AUTHOR::

    Thierry Dumont (2016, 2020).

    EXAMPLES::

    sage: R = PolynomialRing(AA,"x")
    sage: n = 4
    sage: x = P.gen()
    sage: c = [(s[0]+1)/2 for s in R(legendre_P(n,x)).roots()]
    sage: A,B = colloc(c,R,"Gauss 4")
    """
    
    Pb= P.base()
    x = P.gen()
    n = len(c)
    #
    for s in c:
        if s <0 or s >1:
            raise CollocPointNotGood(
                "colloc: collocation point ",s," is out of [0,1]")
    pols=[]
    for  i in range(0,n):
        ploc=P(1)
        for p in range(0,n):
            if i!=p:
                ploc *= P((x-c[p])/(Pb(c[i])-Pb(c[p])))
        pols.append(ploc)

    prims = [p.integral(x) for p in pols]
    prims0 = [p(x = 0) for p in prims]
    A = matrix(Pb,[ [prims[j](x = c[i]) - prims0[j] for j in range(0,n)] \
                  for i in range(0,n)])
    
    B = [prims[j](x = 1) - prims0[j] for j in range(0,n)]
    # exactify to improve lisibility, if possible!
    if R is AA or R is QQbar:
        for i in range(0,n):
            for j in range(0,n):
                A[i,j].exactify()
            B[i].exactify()
    B=vector(B)

    def constructor(self):
        # this will be the contructor of the class returned below.
        self.Title=title
        self.A = A
        self.B = vector(B)
        RungeKutta.__init__(self,A,B,self.Title)
        
    return type("Colloc"+str(len(c)),(RungeKutta,),{
        "__init__": constructor,
        })

