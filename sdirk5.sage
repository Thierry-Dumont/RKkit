class SDIRK5(SageObject):
    title="SDIRKmethod, order 5"
    # See Hairer and Wanner II, page 100.
    A=matrix(AA,5,5)
    gamm=1/4
    for i in range(5):
        A[i,i]=gamm
    A[1,0]=1/2
    A[2,0]=17/50
    A[3,0]=371/1360
    A[4,0]=25/24
    #
    A[2,1]=-1/25
    A[3,1]=-137/2720
    A[4,1]=-49/48
    #
    A[3,2]=15/544
    A[4,2]=125/16
    #
    A[4,3]=-85/12
    #
    B=vector([25/24,-49/48,125/16,-85/12,1/4])
