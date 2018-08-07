# -*- coding: utf-8 -*-
r"""
Explore properties of Runge-Kutta methods.

AUTHOR:
 
 - Thierry Dumont (2016)

"""

from __future__ import absolute_import
from __future__ import print_function

from sage.structure.sage_object import SageObject
from sage.structure.element import generic_power
from sage.arith.power import generic_power
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.matrix.constructor import matrix,identity_matrix
from sage.rings.all import (AA,QQbar)
from sage.symbolic.ring import SR
from sage.symbolic.constants import I
from sage.functions.log import exp
from sage.functions.other import conjugate
from sage.functions.generalized import sign
from sage.calculus.functional import derivative
from sage.misc.cachefunc import cached_function
from sage.rings.infinity import minus_infinity
from sage.misc.lazy_attribute import *

from RKTrees import *
#
from RKExceptions import *
from RKPolutilities import *
#
class  RKformula(SageObject):
    r"""
    Store coefficients of a Runge-Kutta method, compute and store its
    properties.

    EXAMPLES::

    sage: F = RKformula(A,B)
    
    """
    def __init__(self,A,B,C = []):
        """
        Initilalize ``self``.

        EXAMPLES::
        
        sage: A = matrix(AA,[[5/12,-1/12],[3/4,1/4]])
        sage: B = [3/4,1/4]
        sage: F = RKformula(A,B)
        """
        # All elements of the formula must live in an exact ring:
        assert A.base_ring().is_exact(), \
            "Base ring of A must be exact (like QQ, AA, QQbar,...)"
        assert prod([b.parent().is_exact() for b in B]), \
            "Base ring of B elements must be exact (like QQ, AA, QQbar,...)"
        assert prod([c.parent().is_exact() for c in C]), \
            "Base ring of C elements must be exact (like QQ, AA, QQbar,...)"
        self.s = A.dimensions()[0]
        # test A,B, C dimensions compatibility:
        if  self.s !=  A.dimensions()[1] or self.s != len(B):
            raise DimensionsAreIncompatible(A,B,C)
        if C != [] and len(C) != self.s:
            raise DimensionsAreIncompatible(A,B,C)
        # force coefficients to live in AA:
        self.A = A.change_ring(AA)
        self.B = [ AA(b) for b in B]
        self.C = [ AA(c) for c in C]
        #self.D = A.parent().base()
        self.D = AA
        #self.R = PolynomialRing(self.D, 'z')
        self.R =  PolynomialRing(AA, 'z')
        self.s = A.dimensions()[1]

    @lazy_attribute    
    def n_stages(self):
        """
        Return number of stages of the method.
        """
        return self.s
    @lazy_attribute
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
        D = II-z*self.A
        N = D+z*K
        sf = N.determinant()/D.determinant()
        return sf
    @lazy_attribute    
    def A_is_invertible(self):
        """
        Test if the matrix A part of the Butcher array is invertible.

        EXAMPLES::
        
        sage: F.A_is_invertible()
        
        True
        """
        return self.A.is_invertible()

    @lazy_attribute 
    def is_explicit(self):
        """
        Test if the method is explicit.

        EXAMPLES::

        sage: F.is_explicit()

        False

        """
        ret = self.stability_function.denominator().degree()==0
        return ret
    @lazy_attribute   
    def poles_of_stability_function(self):
        """
        Compute the poles of the stability function, if any.
        """
        if self.is_explicit:
            return []
        else:
            RD = self.stability_function.denominator()
            Poles,ok,nr = roots_checked(RD,QQbar)
            if not ok:
                raise RootsException(nr,RD)
            else:
                return Poles
    @lazy_attribute        
    def real_part_of_poles_all_positive(self):
        """
        Documentation is in the name of this method!

        Returns: (all poles have >=0 real part?) and number of poles==0.
        """

        Poles = self.poles_of_stability_function
        llp = len([s for s in Poles if s[0].real()<0 ])
        llzero = len([s for s in Poles if s[0].real()==0 ])
        return llp==0,llzero
    @lazy_attribute
    def order_of_stability_function(self):
        """
        Documentation is in the name of this method!
        """
        order = 0
        z = self.R.gen()
        while derivative(self.stability_function,z,order)(z = 0)==1: \
              order+= 1
        order-= 1
        return order
    @lazy_attribute
    def module_of_stability_function_squared(self):
        """
        Documentation is in the name of this method!
        """
        return self.stability_function*conjugate(self.stability_function)
    @lazy_attribute
    def stability_function_on_im_axis(self):
        """
        Trace of the stability function on the imaginary axis.
        """
        R = self.stability_function
        P = AA['x']
        x = P.gen()
        RIaxe = R(z = QQbar(I)*x)
        return RIaxe
    @lazy_attribute
    def squared_module_of_stability_function_on_Im(self):
        """
        Square of the module of the trace of the stability function on
        the imaginary axis.
        """
        from RKPolutilities import realpart,impart,conj
        RIaxe = self.stability_function_on_im_axis
        RIaxeN = RIaxe.numerator()
        RIaxeD = RIaxe.denominator()
        m2N = RIaxeN*conj(RIaxeD) 
        m2NR = realpart(m2N)
        m2NI = impart(m2N)
        m2n = generic_power(m2NR,2)+generic_power(m2NI,2)
        m2d = RIaxeD*conj(RIaxeD)
        m2 = m2n/generic_power(m2d,2)
        return m2
    @lazy_attribute
    def is_module_of_stability_function_constant_on_Im(self):
        """
        Documentation is in the name of this method!
        """
        m2 = self.squared_module_of_stability_function_on_Im
        x = m2.parent().gen()
        q = derivative(m2,x)
        result =  q==0
        return result
    @lazy_attribute
    def is_module_of_stability_function_less_than_1(self):
        """
        Is the module of the trace of the stability funtion on the
        imaginary axis <1 ?
        """
        is_const = self.is_module_of_stability_function_constant_on_Im
        if is_const:
            lt1 = True
        else:
            m2 = self.squared_module_of_stability_function_on_Im
            l1,ok,n = roots_checked(m2.numerator()-m2.denominator(),QQbar)
            if not ok:
                raise RootsException(n,m2.numerator()-m2.denominator())
            else:
                l1m1 = [s[0] for s in l1 if s[0].imag()==0 ]
                l1m1realp = [s.real() for s in l1m1 if s >= 0]
                l1r = len(l1m1realp)
                if l1r==1:
                    # module value is 0, only in x == 0.
                    lt1 =  m2(z=1) <1
                else:
                    points = [(l1m1realp[i]+l1m1realp[i+1])/2 \
                              for i in range(0,l1r-1)]
                    points.append(l1m1realp[l1r-1]*3/2)
                    signes = [sign(abs(m2(x = s).real())-1) for s in points]
                    lt1 = len([x for x in signes if x >=0]) == 0
            return lt1
    @lazy_attribute    
    def is_A_stable(self):
        """
        Do we have an A_stable method ?

        See: HW II , second edition, page 43.
        """
        ret = (self.is_module_of_stability_function_constant_on_Im or \
               self.is_module_of_stability_function_less_than_1) and  \
               self.real_part_of_poles_all_positive[0]
        return ret
    @lazy_attribute
    def is_stiffly_accurate(self):
        """
        Do we have a stiffly accurate method ?
        """
        # Proposition 3.8 H-W. TII, page 45.
        A = self.A
        B = self.B
        s1 = self.s-1
        Ls =  self.is_A_stable and \
              (len([1 for j in range(0,self.s) if A[s1,j] != B[j]]) == 0 or \
               len([1 for i in range(0,self.s) if A[i,0] != B[0]])  == 0)
        return Ls
    @lazy_attribute
    def is_L_stable(self):
        """
        Do we have a L_stable method ?
        """
        if not self.is_A_stable:
            ok= False
        else:
            R=self.stability_function
            if R.denominator().degree() > R.numerator().degree():
                ok = True
            elif  R.denominator().degree() <= R.numerator().degree():
                ok = False
            return ok
    @lazy_attribute    
    def is_algebraically_stable(self):
        """
        Is the method algebraically stable ?
        """
        if self.is_explicit or len([s for s in self.B if s<0])!=0:
            As = False
        else:
            As=True
            B = self.B
            A = self.A
            s = self.s
            M = matrix(QQbar,s,s)
            for i in range(0,s):
                for j in range(0,s):
                    if B[i]*A[i,j]+B[j]*A[j,i]-B[i]*B[j]<0:
                        As=False
                        if not As:
                            break
                if not As: break
            return As
    @lazy_attribute    
    def conserve_quadratic_invariants(self):
        """
        Documentation is in the name of this method!
        """
        if self.is_explicit:
            ret =  False
        else:
            B = self.B
            A = self.A
            s = self.s
            M = matrix(QQbar,s,s)
            for i in range(0,s):
                for j in range(0,s):
                    M[i,j] = B[i]*A[i,j]+B[j]*A[j,i]-B[i]*B[j]
                ret = M.is_zero()
            return ret
    def check_order_using_rooted_trees(self,order):
        """
        Check rooted tree at order 'order'.
        """
        self.RTrees = RKTrees(order)
        t = False
        for i in range(1,order+1):
            t = self.RTrees.check_order(self.A,self.B,i)
            if not t:
                break
        return t
    @lazy_attribute
    def stability_on_real_negative_axis(self):
        r"""
        In the case where the method is not A-stable, find the limit
        of stability on the real negative axis.
        """
        if self.is_A_stable:
            ret=minus_infinity
        else:
            p=generic_power(self.stability_function,2)-1
            r=[s[0] for s in sorted(p.numerator().roots(),reverse=True)
               if s[0]<0]
            if len(r)==0:
                ret= minus_infinity
            else:
                ret=r[0]
        return ret
    @lazy_attribute
    def order_using_rooted_trees(self):
        """
        Compute the order of the method using rooted trees.
        """
        o = 1
        while self.check_order_using_rooted_trees(o):
            o+= 1
        return o-1
    @lazy_attribute
    def star_function(self,x,y):
        """
        Compute the star function. This is for drawing the "star" associated
        to the formula.
        """
        Rs = self.stability_function
        s = Rs(x+QQbar(I)*y)/exp(x+I*y)
        star=s*conjugate(s)
        self.properties["star_function"]=star
        return star
    def compute_all_properties(self):
        """
        Compute all possible properties of the formula.
        """
        self.A_is_invertible
        self.is_explicit
        self.stability_function
        self.poles_of_stability_function
        self.real_part_of_poles_all_positive
        self.order_of_stability_function
        self.stability_function_on_im_axis
        self.squared_module_of_stability_function_on_Im
        self.is_module_of_stability_function_constant_on_Im
        self.is_module_of_stability_function_less_than_1
        self.is_A_stable
        self.is_stiffly_accurate
        self.is_L_stable
        self.is_algebraically_stable
        self.conserve_quadratic_invariants
        self.stability_on_real_negative_axis
        self.stability_on_real_negative_axis
        self.order_using_rooted_trees
        #x,y=SR.var("x,y")
        #self.star_function(x,y)        
    def print_all_known_properties(self):
        donot=["A","B","C","D","R","s","RTrees"]
        D=self.__dict__
        for key in D:
            if key not in donot:
                print("-> ",key," :\n",D[key],"\n")
