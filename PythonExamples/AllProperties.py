#!/usr/bin/env python
# coding: utf-8

# ## Make all implemented tests on a Runge-Kutta formula ##

# _RKformula_ has a method  _compute_ _all_ _properties()_.
# 
# Let us, for axample compute "all" the properties of a Gauss formula.

try:
    # In passagemath:
    import passagemath_flint
except ImportError:
    #in sage (not passagemath/sage)
    import sage.all
#
from sage.functions.orthogonal_polys import legendre_P
from sage.rings.qqbar import *
#

from RKkit.rkkit import *
#
R = AA["x"]
x=R.gen()

n=3
c = [(s[0]+1)/2 for s in R(legendre_P(n,x)).roots()]
R=RKcolloc.colloc(c,R,"colloc-"+str(n))
F=RKformula(R())


# The following computation can take time (if the order is high; n=3 is relatively high!)).


F.compute_all_properties()


F.print_all_known_properties()
