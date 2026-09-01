#!/usr/bin/env python
# coding: utf-8

# ### Gaussian formulae ###
#
import passagemath_flint
from sage.rings.qqbar import *
from sage.functions.orthogonal_polys import legendre_P
from RKkit  import *
# We must compute in exact numbers, generally algebraic numbers (AA if real, QQbar else).

R = AA["x"]#PolynomialRing(AA,"x")
x=R.gen()


# Build Gaussian RK methods, find stability properties, compute order using rooted trees, draw stability function and order star.
# 
# *Z!* Computing order with rooted trees is very expensive when n grows!

# In[3]:


for n in range(2,4):
    c = [(s[0]+1)/2 for s in R(legendre_P(n,x)).roots()] #roots of shifted polynomials (collocation points)
    GenRK = RKcolloc.colloc(c,R,"Gauss-"+str(n))# coefficients of the RK formula

    #_RKcolloc.colloc as returned a class (a type). We must instantiate it:_

    G=GenRK()

    F=RKformula(G)
    #start some tests
    print("\nNom :",G.Title)
    print("A-stable: ",F.is_A_stable())
    print("L-stable: ",F.is_L_stable())
    print("order: ",F.order())
    print("stability functions: ",F.stability_function())
    print("stiffly accurate? ",F.is_stiffly_accurate())
    print("algebraically stable?",F.is_algebraically_stable())
    print("conserve quadratic invariants?",F.conserve_quadratic_invariants())
    print("symplectic ?",F.is_Symplectic())
    p=RKplot(F,fill=True,ncurves=2,Enlarge=1)
    p.show()
    q=RKplot(F,fill=True,ncurves=2,type="star",Enlarge=1)
    q.show()

