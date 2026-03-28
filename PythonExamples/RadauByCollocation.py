#!/usr/bin/env python
# 
try:
    # In passagemath:
    import passagemath_flint
except ImportError:
    #in sage (not passagemath/sage)
    import sage.all
#
from RKkit.rkkit import *
from RKkit.methods.formulas import *
from RKkit.rkkit.RKcolloc import *
from sage.rings.qqbar import *
from sage.calculus.functional import diff

R = AA["x"]
x=R.gen()


# 
n=3

c=[s[0] for s in diff(x**(n-1)*(x-1)**(n),x,n-1 ).roots()]
print(c)


# The RK. formula:

GenRK = colloc(c,R,"Radau-"+str(n))
#GenRK is a class, we must instantiate it: GenRK() is the instantiation.
F=RKformula(GenRK())

print(F.A)

print(F.B)

F.compute_all_properties()

F.print_all_known_properties()

R=Radau5()

G=RKformula(Radau5())

G.A


G.compute_all_properties()

print(F.stability_function()==G.stability_function())

G.print_all_known_properties()


G.stability_function()

F.stability_function()

F.B

G.B


# _Very strange!_ We did not find the same A and B as in HW book, but we found the same stability function.

# ### Now, order 3 radau method: ###

R = PolynomialRing(AA,"x")
n=2
c=[s[0] for s in diff(x**(n-1)*(x-1)**(n),x,n-1 ).roots()]
print(c)

GenRK = colloc(c,R,"Radau-3")
F=RKformula(GenRK())

F.compute_all_properties()


F.order()

R=RKformula(Radau2a())


print(R.A==F.A and R.B==F.B)





