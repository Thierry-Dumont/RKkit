from sage.all import *
from rkkit.RKRungeKutta import *
class Lobatto4(RungeKutta):
    def __init__(self):
        title="Lobatto method, order 4"
        # Lobatto order 4.See Hairer Nörsett, Wanner TI, page 211.
        A=matrix(AA,[[0,0,0],[1/AA(4),1/AA(4),0],[0,1,0]])
        B=vector(AA,[1/AA(6),AA(2)/AA(3),1/AA(6)])
        super().__init__(A,B,title)
