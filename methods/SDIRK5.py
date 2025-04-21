from sage.all import *
from rkkit.RKRungeKutta import *
class SDIRK5(RungeKutta):
    def __init__(self):
        title="SDIRKmethod, order 5"
        # See Hairer and Wanner II, page 100.
        A=matrix(AA,5,5)
        gamm=1/AA(4)
        for i in range(5):
            A[i,i]=gamm
        A[1,0]=1/AA(2)
        A[2,0]=17/AA(50)
        A[3,0]=371/AA(1360)
        A[4,0]=25/AA(24)
        #
        A[2,1]=-1/AA(25)
        A[3,1]=-137/AA(2720)
        
        A[4,1]=-49/AA(48)
        #
        A[3,2]=15/AA(544)
        A[4,2]=125/AA(16)
        #
        A[4,3]=-85/AA(12)
        #
        B=vector([25/AA(24),-49/AA(48),125/AA(16),-85/AA(12),1/AA(4)])
        super().__init__(A,B,title)
