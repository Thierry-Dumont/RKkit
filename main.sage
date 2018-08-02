from RKformula import *
from RKplot import *
from RKcolloc import *
load("radau5.sage")
F=RKformula(A,B)
R = PolynomialRing(AA,"x")
n =2
c = [(s[0]+1)/2 for s in R(legendre_P(n,x)).roots()]
A,B = colloc(c,R)
F=RKformula(A,B)