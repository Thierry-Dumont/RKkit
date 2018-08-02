class Radau5:
    s=3#nb stages.
    title="Radau5 (IIA) method"
    comment= title
    A=matrix(AA,s,s)
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
    B=[1/9,(16+s6)/36,(16-s6)/36]
