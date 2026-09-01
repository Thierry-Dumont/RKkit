#!/usr/bin/env python
# coding: utf-8

# ## Test: can we compute all  properties for all coded formulas? ##
import passagemath_flint
from sage.rings.qqbar import *
import sage.rings.qqbar
#
from RKkit import *
from RKkit.formulas import *
from RKkit.RKcolloc import *

cformulas=[Radau5(),Gauss4(),Radau2a(),RK4(),SDIRK5(),Lobatto4()]

# ,SDIRK3()


for formula in cformulas:
    print(formula.Title)
    F=RKformula(formula)
    F.compute_all_properties()
    print("* ",formula.Title,': ok.\n')
    ccc=input("type any character to continue:")

# Try also some Gauss formulas.



R = PolynomialRing(AA,"x")
for n in[1,2,3]:
    x = R.gen()
    c = [(s[0]+1)/2 for s in R(legendre_P(n,x)).roots()]
    GenRK = colloc(c,R,"Gauss-"+str(n))
    G=GenRK()
    F=RKformula(G)
    F.compute_all_properties()
    print("* ",G.Title,': ok.\n')


# ##### end #####






