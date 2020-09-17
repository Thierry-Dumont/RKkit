# -*- coding: utf-8 -*-
"""
Exceptions for the package rkkit.
"""
from sage.all import *
#
class DimensionsAreIncompatible(SageObject,Exception):
    r"""
    Exception raised when the sizes of A , B  and C
    in a Runge-Kutta are not equal. 
    """
    def __init__(self,A,B,C):
        self.An = A.dimensions()[0]
        self.Bn = len(B)
        self.Cn = len(C)
    def __str__(self):
        ret = "Butcher Array: dimension of A and B (or C) are not equal: " \
            +str(self.An)+" , "+str(self.Bn)
        if self.Cn != 0: ret+= " , "+str(self.Cn)
class RootsException(Exception):
    r"""
    Exception raised when the roots of a given polynomial
    could not be *all* computed (with their multiplicity).
    (We generally compute in QQbar, and we can meet Evariste Gallois). 
    """
    def __init__(self,ncomp,pol):
        self.ncomp = ncomp
        self.degree = pol.degree
        self.pol = pol
    def __str__(self):
        return "for " + str(self.pol) + " (degree= " + \
            str(self.degree) + " ), only " \
            +str(self.ncomp)+ " where computed."
class MatrixIsSingular(Exception):
    """
    Exception raised when a matrix is singular/
    """
    def __init__(self,t):
        self.t = t
    def __str__(self):
        return "for this method, "+self.t+ " is a singular matrix"
class GraphicProblem(Exception):
    """
    Exception raised when an error happens in a graphic.
    """
    def __init__(self,t):
        self.t = t
    def __str__(self):
        return "Graphic problem: "+self.t

