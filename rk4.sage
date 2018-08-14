class RK4(SageObject):
    title="Classical Runge-Kutta 4 explicit method"
    A=matrix(AA,[[0,0,0,0],[1/2,0,0,0],[0,1/2,0,0],[0,0,1,0]])
    B=vector([1/6,2/6,2/6,1/6])
