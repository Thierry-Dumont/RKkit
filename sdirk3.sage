from sage.all import *
class SDIRK3(SageObject):
    g=AA((2+sqrt(3))/6)
    A=matrix(AA,[[g,0],[1-2*g,g]])
    B=vector([1/2,1/2])
    title= "Sdirk3 method"
    
