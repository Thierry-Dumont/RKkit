from sage.all import *
from rkkit.RKRungeKutta import *
class Radau2a(RungeKutta):
    def __init__(self):
        title="Radau 2a method"
        A=matrix(AA,[[5/AA(12),-1/AA(12)],[3/AA(4),1/AA(4)]])
        ##
        B=vector([3/AA(4),1/AA(4)])
        super().__init__(A,B,title)
