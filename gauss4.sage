from sage.all import *
from rkkit.RKRungeKutta import *
class Gauss4(RungeKutta):
    def __init__(self):
        title="Gauss method, order 2."
        A=matrix(AA,[[1/4,1/4-sqrt(3)/6],[1/4+sqrt(3)/6,1/4]])
        ##
        B=vector([1/2,1/2])
        super().__init__(A,B,title)
