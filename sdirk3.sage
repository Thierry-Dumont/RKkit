from sage.all import *
from rkkit.RKRungeKutta import *
class SDIRK3(RungeKutta):
    def __init__(self):
        g=AA((2+sqrt(3))/6)
        A=matrix(AA,[[g,0],[1-2*g,g]])
        B=vector([1/2,1/2])
        title= "Sdirk3 method"
        super().__init__(A,B,title)
    
