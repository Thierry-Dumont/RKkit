from sage.all import *
from .RKExceptions import *
class RungeKutta(SageObject):
    """
    Base class for all Runge-Kutta methods.

    """
    def __init__(self,A,B,Title,C=[]):
        if not A.parent().is_exact():
            raise MustBeExact("RungeKutta: parent of A is not exact")
        if not  B.parent().is_exact():
            raise MustBeExact("RungeKutta: parent of B is not exact")
        if not isinstance(A,sage.structure.element.Matrix):
            raise NotA("RungeKutta: A is not a matrix")
        if not isinstance(B,sage.structure.element.Vector):
            raise NotA("RungeKutta: B is not a vector")
        if  A.dimensions()[0] != A.dimensions()[1]\
            or A.dimensions()[0] != len(B):
            raise  DimensionsAreIncompatible(A,B,C)
        if C != [] and len(C) != A.dimensions()[0]:
            raise DimensionsAreIncompatible(A,B,C)
            
        
        self.A = A
        self.B = B
        self.Title = Title
        self.C = C
    def __str__(self):
        return str(self.A)+"\n"+str(self.B)
