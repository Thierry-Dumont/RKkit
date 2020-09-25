from sage.all import *
from rkkit.RKRungeKutta import *
class Lobatto4(RungeKutta):
    def __init__(self):
        title="Lobatto method, order 4"
        # Lobatto order 4.See Hairer Nörsett, Wanner TI, page 211.
        A=matrix(AA,[[0,0,0],[1/4,1/4,0],[0,1,0]])
        B=vector([1/6,2/3,1/6])
        super().__init__(A,B,title)
