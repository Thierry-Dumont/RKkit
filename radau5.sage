from sage.all import *
from rkkit.RKRungeKutta import *
class Radau5(RungeKutta):
    def __init__(self):
        title="Radau5 method"
        #
        A=matrix(AA,3,3)
        s6=AA(sqrt(6))
        A[0,0]=1/9
        A[0,1]=(-1-s6)/18
        A[0,2]=(-1+s6)/18
        #
        A[1,0]=1/9
        A[1,1]=(88+7*s6)/360
        A[1,2]=(88-43*s6)/360
        #
        A[2,0]=1/9
        A[2,1]=(88+43*s6)/360
        A[2,2]=(88-7*s6)/360
        ##
        B=vector([1/9,(16+s6)/36,(16-s6)/36])
        super().__init__(A,B,title)

