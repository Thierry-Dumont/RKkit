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

class Lobatto4(RungeKutta):
    def __init__(self):
        title="Lobatto method, order 4"
        # Lobatto order 4.See Hairer Nörsett, Wanner TI, page 211.
        A=matrix(AA,[[0,0,0],[1/AA(4),1/AA(4),0],[0,1,0]])
        B=vector(AA,[1/AA(6),AA(2)/AA(3),1/AA(6)])
        super().__init__(A,B,title)

class Radau2a(RungeKutta):
    def __init__(self):
        title="Radau 2a method"
        A=matrix(AA,[[5/AA(12),-1/AA(12)],[3/AA(4),1/AA(4)]])
        ##
        B=vector([3/AA(4),1/AA(4)])
        super().__init__(A,B,title)

class Radau5(RungeKutta):
    def __init__(self):
        title="Radau5 method"
        #
        A=matrix(AA,3,3)
        s6=AA(sqrt(6))
        A[0,0]=1/AA(9)
        A[0,1]=(-1-s6)/18
        A[0,2]=(-1+s6)/18
        #
        A[1,0]=1/AA(9)
        A[1,1]=(88+7*s6)/360
        A[1,2]=(88-43*s6)/360
        #
        A[2,0]=1/AA(9)
        A[2,1]=(88+43*s6)/360
        A[2,2]=(88-7*s6)/360
        ##
        B=vector([1/AA(9),(16+s6)/36,(16-s6)/36])
        super().__init__(A,B,title)


class RK4(RungeKutta):
    def __init__(self):
        title="Classical Runge-Kutta 4 explicit method"
        A=matrix(AA,[[0,0,0,0],[1/AA(2),0,0,0],
                     [0,1/AA(2),0,0],[0,0,1,0]])
        B=vector(AA,[1/AA(6),2/AA(6),2/AA(6),1/AA(6)])
        super().__init__(A,B,title)

class SDIRK3(RungeKutta):
    def __init__(self):
        g=AA((2+sqrt(3))/6)
        A=matrix(AA,[[g,0],[1-2*g,g]])
        B=vector([1/AA(2),1/AA(2)])
        title= "Sdirk3 method"
        super().__init__(A,B,title)
    

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
class veryBad(RungeKutta):
    def __init__(self):
        title="Classical Runge-Kutta 4 explicit method, ill coded"
        # 1/2,2/6 ... will be converted into floats, which are not exact!
        A=matrix(AA,[[0,0,0,0],[1/2,0,0,0],
                     [0,1/2,0,0],[0,0,1,0]])
        B=vector(AA,[1/6,2/6,2/6,1/6])
        super().__init__(A,B,title)
