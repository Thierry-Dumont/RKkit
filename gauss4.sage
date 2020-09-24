from sage.all import *
class Gauss4(SageObject):
    title="Gauss method, order 2."
    A=matrix(AA,[[1/4,1/4-sqrt(3)/6],[1/4+sqrt(3)/6,1/4]])
    ##
    B=vector([1/2,1/2])
