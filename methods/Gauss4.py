from sage.all import *
from rkkit.RKRungeKutta import *
class Gauss4(RungeKutta):
    def __init__(self):
        title="Gauss method, 2 points, order 4."
        A=matrix(AA,[[1/AA(4),1/AA(4)-sqrt(3)/AA(6)],[1/AA(4)+sqrt(3)/AA(6),
                                                      1/AA(4)]])
        ##
        B=vector([1/AA(2),1/AA(2)])
        super().__init__(A,B,title)
