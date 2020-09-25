from sage.all import *
from .RKExceptions import *
#from sage.matrix.matrix import is_Matrix
class RungeKutta(SageObject):
    def __init__(self,A,B,Title,C=[]):
        if not A.parent().is_exact():
            raise Exception("RungeKutta: parent of A is not exact")
        if not  B.parent().is_exact():
            raise Exception("RungeKutta: parent of B is not exact")
        # print(A.parent())
        # if not is_Matrix(A):
        #     print("A=",A)
        #     #raise Exception("RungeKutta: A is not a matrix")
        # if not is_Matrix(B):
        #     print("B=",B)
        #     #raise Exception("RungeKutta: A is not a matrix")
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
