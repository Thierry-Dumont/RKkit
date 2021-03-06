# -*- coding: utf-8 -*-
r"""

Explore properties of Runge-Kutta methods.

AUTHOR:
 
 - Thierry Dumont (2016, 2018. 2020: Python 3, _persistance decorator)
   Institut C. Jordan, Lyon, France.
"""

# ****************************************************************************
#       Copyright (C) 2018 Thierry Dumont tdumont@math.univ-lyon1.fr
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.all import *
from sage.rings.infinity import minus_infinity
#
from .RKTrees import *
#
from .RKExceptions import *
from .RKPolutilities import *
#
import functools
#
class  RKformula(SageObject):
    r"""
    Store coefficients of a Runge-Kutta method, compute some of its
    properties.

    EXAMPLES::

    sage: F = RKformula(A,B)
    
    """
    def __init__(self,F):
        """
        Initilalize ``self``. F is a Runge-Kutta class.

        EXAMPLES::
        
        sage: R = RK4()
        sage: F = RKformula(R)
        """
       
        # force coefficients to live in AA:
        self.A = F.A.change_ring(AA)
        self.B = vector([ AA(b) for b in F.B])
        self.C = vector([ AA(c) for c in F.C])
        self.D = AA
        self.R =  PolynomialRing(AA, 'z')
        self.s = self.A.dimensions()[1]
        # as computing properties can be slow, we will cache them here as
        # soon as they are computed:
        self.known_properties={}

           
    def _persistance(foo):
        """
        Decorator: caches results of "foo" in self.known_properties
        """
        @functools.wraps(foo)
        def magic( self, *args, **kwargs ):
            if foo.__name__ in self.known_properties:
                return self.known_properties[foo.__name__]
            else:
                x=foo( self,  *args, **kwargs )
                self.known_properties[foo.__name__]=x
                return x
        return magic
    
    def _latex_(self):
        r"""
        Return the LaTeX representation of X.
        """
        return latex(self.A)+" "+latex(self.B)

    def n_stages(self):
        """
        Return number of stages.
        """
        return self.s
    
    @_persistance
    def stability_function(self):
        """

        Compute the stability function.

        EXAMPLES::

        sage: F.stability_function()

        (1/3*z + 1)/(1/6*z^2 - 2/3*z + 1)
        """
        z = self.R.gen()
        Rng = self.D
        K = matrix(Rng,[self.B for i in range(0,self.s)])
        II = identity_matrix(Rng,self.s)
        D = II - z*self.A
        N = D + z*K
        return N.determinant()/D.determinant()

    @_persistance    
    def A_is_invertible(self):
        """

        Test if the matrix A part of the Butcher array is invertible.

        EXAMPLES::
        
        sage: F.A_is_invertible()
        
        True
        """
        return self.A.is_invertible()

    @_persistance 
    def is_explicit(self):
        """
        Test if the method is explicit.

        EXAMPLES::

        sage: F.is_explicit()

        False

        """
        return self.stability_function().denominator().degree()==0

    @_persistance   
    def poles_of_stability_function(self):
        """
        Compute the poles of the stability function, if any.
        """
        if self.is_explicit():
            return []
        else:
            RD = self.stability_function().denominator()
            Poles,ok,nr = roots_checked(RD,QQbar)
            if not ok:
                raise RootsException(nr,RD)
            else:
                return Poles
            
    @_persistance        
    def real_part_of_poles_all_positive(self):
        """
        Documentation is in the name of this method!

        Returns: (all poles have >=0 real part?) and number of poles==0.
        """
        Poles = self.poles_of_stability_function()
        llp = all(s[0].real()>=0 for s in Poles)
        llzero = len([s for s in Poles if s[0].real()==0 ])
        return llp,llzero
    
    @_persistance
    def order_of_stability_function(self):
        """
        Documentation is in the name of this method!
        """
        order = 0
        z = self.R.gen()
        while derivative(self.stability_function(),z,order)(z = 0)==1: \
              order+= 1
        return order -1
    
    @_persistance
    def module_of_stability_function_squared(self):
        """
        Documentation is in the name of this method!
        """
        return self.stability_function()*conjugate(self.stability_function())
    
    @_persistance
    def stability_function_on_im_axis(self):
        """
        Trace of the stability function on the imaginary axis.
        """
        R = self.stability_function()
        P = AA['x']
        x = P.gen()
        return R(z = QQbar(I)*x)

    @_persistance
    def squared_module_of_stability_function_on_Im(self):
        """
        Square of the module of the trace of the stability function on
        the imaginary axis.
        """
        RIaxe = self.stability_function_on_im_axis()
        RIaxeN = RIaxe.numerator()
        RIaxeD = RIaxe.denominator()
        m2N = RIaxeN*conj(RIaxeD) 
        m2NR = realpart(m2N)
        m2NI = impart(m2N)
        m2n = generic_power(m2NR,2)+generic_power(m2NI,2)
        m2d = RIaxeD*conj(RIaxeD)
        return m2n/generic_power(m2d,2)
    
    @_persistance
    def is_module_of_stability_function_constant_on_Im(self):
        """
        Documentation is in the name of this method!
        """
        m2 = self.squared_module_of_stability_function_on_Im()
        x = m2.parent().gen()
        return derivative(m2,x) == 0

    @_persistance
    def is_module_of_stability_function_less_than_1(self):
        """
        Is the module of the trace of the stability funtion on the
        imaginary axis <1 ?
        """
        is_const = self.is_module_of_stability_function_constant_on_Im()
        if is_const:
            return True
        else:
            m2 = self.squared_module_of_stability_function_on_Im()
            l1,ok,n = roots_checked(m2.numerator()-m2.denominator(),QQbar)
            if not ok:
                raise RootsException(n,m2.numerator()-m2.denominator())
            else:
                l1m1 = [s[0] for s in l1 if s[0].imag()==0 ]
                l1m1realp = [s.real() for s in l1m1 if s >= 0]
                l1r = len(l1m1realp)
                if l1r==1:
                    # module value is 0, only in x == 0.
                    return  m2(z=1) <1
                else:
                    points = [(l1m1realp[i]+l1m1realp[i+1])/2 \
                              for i in range(0,l1r-1)]
                    points.append(l1m1realp[l1r-1]*3/2)
                    signes = [sign(abs(m2(x = s).real())-1) for s in points]
                    return len([x for x in signes if x >=0]) == 0
    
    @_persistance    
    def is_A_stable(self):
        """
        Do we have an A_stable method ?

        See: HW II , second edition, page 43.
        """
        return  (self.is_module_of_stability_function_constant_on_Im() or \
                 self.is_module_of_stability_function_less_than_1()) and  \
                 self.real_part_of_poles_all_positive()[0]

    @_persistance
    def is_stiffly_accurate(self):
        """
        Do we have a stiffly accurate method ?
        """
        # Proposition 3.8 in H-W. TII, page 45.
        A = self.A
        B = self.B
        s1 = self.s-1
        return self.is_A_stable() and \
              (len([1 for j in range(0,self.s) if A[s1,j] != B[j]]) == 0 or \
               len([1 for i in range(0,self.s) if A[i,0] != B[0]])  == 0)

    @_persistance
    def R_infinite(self):
        r"""
        R(\infty). See Hairer-Wanner t.II pages 45 and 375.
        """
        if self.A_is_invertible():
            AI = self.A.inverse()
            One = vector([AA(1) for i in range(0,self.s)])
            ret = AA(1)-self.B.dot_product(AI*One)
            ret.exactify()
            return ret
        else:
            raise MatrixIsSingular("A")
    
    @_persistance
    def is_L_stable(self):
        """
        Do we have a L_stable method ?
        """
        if not self.is_A_stable():
            return False
        else:
            R=self.stability_function()
            return R.denominator().degree() > R.numerator().degree() or \
                self.R_infinite().is_zero()
        

    @_persistance
    def M_matrix(self):
        """

        Compute the matrix M used by is_algebraically_stable(self) and
        conserve_quadratic_invariants(self).

        """
        B = self.B
        A = self.A
        M = matrix(QQbar,self.s,self.s)
        for i in range(0,self.s):
            for j in range(0,self.s):
                M[i,j] = B[i]*A[i,j]+B[j]*A[j,i]-B[i]*B[j]
        return M
    
    @_persistance    
    def is_algebraically_stable(self):
        """

        Is the method algebraically stable ?

        """
        
        if self.is_explicit() or any(s<0 for s in self.B):
            return False
        else:
            M = self.M_matrix()
            return all(s >=0 for s in M.list())
        
    @_persistance
    def is_Symmetric(self):
        """
        All is in the title.
        """
        P = matrix(AA,self.s,self.s)
        for i in range(0,self.s):
            P[i,self.s-i-1]=1
        PA =P *self.A+ self.A*P
        return all(self.B == Row for Row in PA.rows()) and \
            self.B == P*self.B

    @_persistance    
    def conserve_quadratic_invariants(self):
        """
        Documentation is in the name of this method!
        """
        return not self.is_explicit() and self.M_matrix().is_zero()
    
    @_persistance
    def is_Symplectic(self):
        """

        Test for simplexity. We can conclude (positively) only if
        the matrix N (see above) is zero.

        """
        if self.is_explicit():
            return False
        elif self.M_matrix().is_zero():
            return True
        else:
            return "Unknown"
        
    def check_order_using_rooted_trees(self,order):
        """

        Check rooted tree at order 'order'.

        """
        self.RTrees = RKTrees()
        t = False
        for i in range(1,order+1):
            t = self.RTrees.check_order(self.A,self.B,i)
            if not t:
                break
        return t #True iff formula is at least of order "order".
    
    @_persistance
    def make_order_equations(self,order):
        self.RTrees = RKTrees()
        s=[]
        for i in range(1,order+1):
            s+=self.RTrees.make_order_equations(self.A,self.B,i)
        return s
    @_persistance
    def stability_on_real_negative_axis(self):
        """

        In the case where the method is not A-stable, find the limit
        of stability on the real negative axis.

        """
        if self.is_A_stable():
            return minus_infinity
        else:
            p=generic_power(self.stability_function(),2)-1
            r=[s[0] for s in sorted(p.numerator().roots(),reverse=True)
               if s[0]<0]
            if len(r)==0:
                return minus_infinity
            else:
                return r[0]
 
    @_persistance
    def order(self):
        """

        Compute order of the method using rooted trees.

        """
        o = 0
        while self.check_order_using_rooted_trees(o+1):
            o+= 1
        return o
    
    @_persistance
    def order_star_function(self):
        """
        Compute the order star function. 
        This is for drawing the "star" associated to the formula.
        """
        x=SR.var("x")
        y=SR.var("y")
        Rs = self.stability_function()
        RRs=Rs.numerator().change_ring(RDF)/Rs.denominator().change_ring(RDF)
        s=RRs(x+SR(I)*y)/exp(x+SR(I)*y)
        return s.abs()

    def compute_all_properties(self):
        """

        Compute all possible properties of the formula.

        """
        self.A_is_invertible()
        self.is_explicit()
        self.stability_function()
        self.poles_of_stability_function()
        self.real_part_of_poles_all_positive()
        self.order_of_stability_function()
        self.stability_function_on_im_axis()
        self.squared_module_of_stability_function_on_Im()
        self.is_module_of_stability_function_constant_on_Im()
        self.is_module_of_stability_function_less_than_1()
        self.is_A_stable()
        self.is_stiffly_accurate()
        self.is_L_stable()
        self.is_algebraically_stable()
        self.is_Symmetric()
        self.is_Symplectic()
        self.conserve_quadratic_invariants()
        self.stability_on_real_negative_axis()
        self.stability_on_real_negative_axis()
        self.order()
        self.order_star_function()
        
    def print_all_known_properties(self):
        """

        Print all already computed properties.

        """
        donotprint=["A","B","C","D","R","s","RTrees","M_matrix"]
        D=self.known_properties
        for key in D:
            if key not in donotprint:
                print("-> ",key," : ",D[key],"\n")
