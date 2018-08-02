# -*- coding: utf-8 -*-
from __future__ import absolute_import 
from sage.structure.sage_object import SageObject
#from sage.structure.element import generic_power
from sage.arith.power  import generic_power
from sage.rings.all import (QQ,AA,QQbar)
from sage.symbolic.ring import SR
from sage.symbolic.constants import I
def conj(P):
    # conjugue d'un polynome:
    Pc = P.coefficients()
    x = P.parent().gen()
    return sum([Pc[p]*generic_power(x,p) for p in range(0,len(Pc))\
                if Pc[p].imag() == 0]) - \
                   sum([Pc[p]*generic_power(x,p) for p in range(0,len(Pc)) \
                        if Pc[p].imag() != 0])
def impart(P):
    # partie imaginaire:
    Pc = P.coefficients()
    x = P.parent().gen() 
    Im = QQbar(I)
    return sum([-Im*Pc[p]*generic_power(x,p) \
                for p in range(0,len(Pc))  if Pc[p].imag() != 0])
def realpart(P):
    # partie reelle:
    x = P.parent().gen()
    Pc = P.coefficients()
    return sum([Pc[p]*generic_power(x,p) for p in range(0,len(Pc)) \
                if Pc[p].imag() == 0])
def roots_checked(pol,Q):
    q = pol.change_ring(Q)
    rac = q.roots()
    n = sum([s[1] for s in rac])
    return rac,pol.degree()== n, n
