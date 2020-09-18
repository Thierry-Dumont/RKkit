# -*- coding: utf-8 -*-
"""
Implement somme non implemented computations (in Sgae) for polynomials.
"""
from sage.all import *
from sage.arith.power  import generic_power
#
def conj(P):
    """
    Conjugate of a polynomial P.
    """
    Pc = P.coefficients()
    x = P.parent().gen()
    return sum([Pc[p] * generic_power(x,p) for p in range(0,len(Pc))\
                if Pc[p].imag() == 0]) - \
                   sum([Pc[p] * generic_power(x,p) for p in range(0,len(Pc)) \
                        if Pc[p].imag() != 0])
def impart(P):
    """
    Imaginary part of a polynomial P.

    (ie: compute sum Im P_i x^i only  for P_i imaginary).

    """
    Pc = P.coefficients()
    x = P.parent().gen() 
    Im = QQbar(I)
    return sum([-Im * Pc[p]*generic_power(x,p) \
                for p in range(0,len(Pc))  if Pc[p].imag() != 0])
def realpart(P):
    """
    Real part of a polynomial P.

    (ie: compute sum P_i x^i  only for P_i real).

    """
    x = P.parent().gen()
    Pc = P.coefficients()
    return sum([Pc[p] * generic_power(x,p) for p in range(0,len(Pc)) \
                if Pc[p].imag() == 0])
def roots_checked(pol,R):
    """
    Check that, computing in the ring R, we can compute n roots (with their
    multiplicity) for a  polynomial pol of degree n.
    """
    q = pol.change_ring(R)
    rac = q.roots()
    n = sum([s[1] for s in rac])
    return rac,pol.degree()== n, n

