#!/usr/bin/env python
# coding: utf-8

# # Playing with Runge-Kutta methods ##

import sage.all
#
from RKkit.rkkit import *
# Let us import somme predefined Runge-Kutta method descriptions:
from RKkit.methods.formulas import *


# Choose a formula:

#RK=Radau5()
RK=Lobatto4()
#RK=SDIRK3()
#RK=Gauss4()
#RK=Radau2a()
#RK=RK4()
#RK=SDIRK5()

print(RK)

# and define the formula:
F=RKformula(RK)


# Ok, now let-us check different properties of the formula:

print(F.is_explicit())

print(F.is_A_stable())

print(F.is_L_stable())

print(F.is_stiffly_accurate())


print(F.is_algebraically_stable())


print(F.conserve_quadratic_invariants())



print(F.stability_function())


print(F.poles_of_stability_function())


# Find the limit of the stability domain on $\mathbb{R}^-$. It can be $-\infty$, for example if the method is A-stable:


print(F.stability_on_real_negative_axis())


print(F.order_of_stability_function())


# The true order (computed using rooted trees):



print(F.order())


# Symplectic ?

print(F.is_Symplectic())


# Now, let us plot the stability domain:

# In[17]:


#@interact
def P(WindowSize=(0.1,10,0.5),Translate=(-50,50,5)):
    RKplot(F,fill=True,ncurves=2,Enlarge=WindowSize,TranslateX=Translate).show()


# We can also plot the orde star (slow):


#@interact
def PStar(WindowSize=(0.01,10,0.5),Translate=(-50,50,5)):
    RKplot(F,fill=True,ncurves=1,type="star",Enlarge=WindowSize,TranslateX=Translate).show()


# In[ ]:




