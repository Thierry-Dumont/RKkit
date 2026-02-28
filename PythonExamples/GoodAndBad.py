#!/usr/bin/env python
# coding: utf-8

# ## This small notebook illustrates how to create a Runge-Kutta method.
# 
# Remember that we must use exact sets of number, namely (real) algebraic nubers (AA or QQbar) or rationals (QQ) when defining Runge-Kutta methods.


from RKkit.rkkit import *
from RKkit.rkkit.RKRungeKutta import *
from RKkit.methods.formulas import *
from sage.rings.qqbar import *

#  A "good" initialization (this the classical explicit RK4 mathod).
#  
#  Here all numbers are real algebraic numbers, or rational numbers, but __not__ floating *point numbers*: we need to make exact computations. 


class GoodRK4(RungeKutta):
    def __init__(self):
        title="Classical Runge-Kutta 4 explicit method"
        A=matrix(AA,[[0,0,0,0],[1/AA(2),0,0,0],
                     [0,1/AA(2),0,0],[0,0,1,0]])
        B=vector(AA,[1/AA(6),2/AA(6),2/AA(6),1/AA(6)])
        super().__init__(A,B,title)



Good=RKformula(GoodRK4())


# So, all went well.
cc=input('"GoodRK4" ok! type any character to continue with the "bad" coded method')

# A "bad" initialization, but which could appear as "good":


class NotSoBadRK4(RungeKutta):
    def __init__(self):
        title="Classical Runge-Kutta 4 explicit method"
        A=matrix(AA,[[0,0,0,0],[1/2,0,0,0],
                     [0,1/2,0,0],[0,0,1,0]])
        B=vector(AA,[1/6,2/6,2/6,1/6])
        super().__init__(A,B,title)




NotSoBad=RKformula(NotSoBadRK4())


# Let us look at the "A" matrix:


NotSoBad.A



NotSoBad.A.parent()


# Why that ? 

# -Here everything is ok, since "NotSoBad" coefficients was transformed by the sage _preparser,_ which avoid creation of floats.
# 
# -But if we import the "veryBad" from methods.formulas, which is a python file (have a look at the bottom of methods/formulas), if is the "same" code as in "NotSoBadRK4", but it is **not** _preparsed_ by sage. 
# 
# Run the next cell and have a look at error message:


v=veryBad()


# The TypeError says that (some) nembers could not be coerced to algebraic numbers (because they where transformed into floats).
# 

# ## For safety, always code your formulas as in "GoodRK4" at the top of this notebook!
# 
# Code 1/AA(2) and not 1/2, for example.
# 




