from sage.all import *
from rkkit.RKRungeKutta import *
class RK4(RungeKutta):
    def __init__(self):
        title="Classical Runge-Kutta 4 explicit method"
        A=matrix(AA,[[0,0,0,0],[1/2,0,0,0],[0,1/2,0,0],[0,0,1,0]])
        B=vector([1/6,2/6,2/6,1/6])
        super().__init__(A,B,title)
