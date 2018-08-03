# -*- coding: utf-8 -*-
r"""
Explore properties of Runge-Kutta methods.

AUTHOR:
 
 - Thierry Dumont (2016)

"""

from __future__ import absolute_import

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
        # 
        self.properties = {}
    def n_stages(self):
        """
        Return number of stages of the method.
        """
        return self.s

    def stability_function(self):
        """
        Compute the stability function.

        EXAMPLES::

        sage: F.stability_function()

        (1/3*z + 1)/(1/6*z^2 - 2/3*z + 1)
        """
        if "Stability_Function" in self.properties:
            return self.properties["Stability_Function"]
        else:
            z = self.R.gen()
            Rng = self.D
            K = matrix(Rng,[self.B for i in range(0,self.s)])
            II = identity_matrix(Rng,self.s)
            D = II-z*self.A
            N = D+z*K
            sf =  N.determinant()/D.determinant()
            self.properties["Stability_Function"] = sf
            return sf
    def A_is_invertible(self):
        """
        Test if the matrix A part of the Butcher array is invertible.

        EXAMPLES::
        
        sage: F.A_is_invertible()
        
        True
        """
        if "A_is_invertible" in self.properties:
             return self.properties["A_is_invertible"]
        else:
            ok = self.A.is_invertible()
            self.properties["A_is_invertible"] = ok
            return ok
    def is_explicit(self):
        """
        Test if the method is explicit.

        EXAMPLES::

        sage: F.is_explicit()

        False

        """
        if "is_explicit" in self.properties:
            return self.properties["is_explicit"]
        else:
            ret = self.stability_function().denominator().degree()==0
            self.properties["is_explicit"] = ret
            return ret
        
    def poles_of_stability_function(self):
        """
        Compute the poles of the stability function, if any.
        """
        if "Poles" in  self.properties:
            return self.properties["Poles"]
        else:
            if self.is_explicit():
                self.properties["Poles"] = []
                return []
            else:
                RD = self.stability_function().denominator()
                Poles,ok,nr = roots_checked(RD,QQbar)
                if not ok:
                    raise RootsException(nr,RD)
                else:
                    self.properties["Poles"] = Poles
                    return Poles
    def real_part_of_poles_all_positive(self):
        """
        Documentation is in the name of this method!
        """
        if "real_part_of_poles_all_positive" in self.properties:
            llp = self.properties["real_part_of_poles_all_positive"][0]
            llzero = self.properties["real_part_of_poles_all_positive"][1]
            return llp,llzero
        else:
            Poles = self.poles_of_stability_function()
            llp = len([s for s in Poles if s[0].real()<0 ])
            llzero = len([s for s in Poles if s[0].real()==0 ])
            self.properties["real_part_of_poles_all_positive"] = (llp==0,llzero)
            return llp==0,llzero
   
    def order_of_stability_function(self):
        """
        Documentation is in the name of this method!
        """
        if "order_of_stability_function" in  self.properties:
            return  self.properties["order_of_stability_function"]
        else:
            order = 0
            z = self.R.gen()
            while derivative(self.stability_function(),z,order)(z = 0)==1: \
                  order+= 1
            order-= 1
            self.properties["order_of_stability_function"] = order
            return order
    def module_of_stability_function_squared(self):
        """
        Documentation is in the name of this method!
        """
        return self.stability_function()*conjugate(self.stability_function())
    def stability_function_on_im_axis(self):
        """
        Trace of the stability function on the imaginary axis.
        """
        if "stability_function_on_im_axis" in self.properties:
            return self.properties["stability_function_on_im_axis"]
        else:
            R = self.stability_function()
            P = AA['x']
            x = P.gen()
            RIaxe = R(z = QQbar(I)*x)
            self.properties["stability_function_on_im_axis"] = RIaxe
            return RIaxe
    def squared_module_of_stability_function_on_Im(self):
        """
        Square of the module of the trace of the stability function on
        the imaginary axis.
        """
        if "squared_module_of_stability_function_on_Im" in self.properties:
            return self.properties["squared_module_of_stability_function_on_Im"]
        else:
            from RKPolutilities import realpart,impart,conj
            RIaxe = self.stability_function_on_im_axis()
            RIaxeN = RIaxe.numerator()
            RIaxeD = RIaxe.denominator()
            m2N = RIaxeN*conj(RIaxeD) 
            m2NR = realpart(m2N)
            m2NI = impart(m2N)
            m2n = generic_power(m2NR,2)+generic_power(m2NI,2)
            m2d = RIaxeD*conj(RIaxeD)
            m2 = m2n/generic_power(m2d,2)
            self.properties["squared_module_of_stability_function_on_Im"] = m2
            return m2
    def is_module_of_stability_function_constant_on_Im(self):
        """
        Documentation is in the name of this method!
        """
        if "is_module_of_stability_function_constant_on_Im" in self.properties:
            return \
                self.properties["is_module_of_stability_function_constant_on_Im"]
        else:
            m2 = self.squared_module_of_stability_function_on_Im()
            x = m2.parent().gen()
            q = derivative(m2,x)
            result =  q==0
            self.properties["is_module_of_stability_function_constant_on_Im"] \
                = result
            return result
    def is_module_of_stability_function_less_than_1(self):
        """
        Is the module of the trace of the stability funtion on the
        imaginary axis <1 ?
        """
        if "is_module_of_stability_function_less_than_1" in self.properties:
            return self.properties["is_module_of_stability_function_less_than_1"]
        else:
            is_const = self.is_module_of_stability_function_constant_on_Im()
            if is_const:
                lt1 = True
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
                        lt1 =  m2(z=1) <1
                    else:
                        points = [(l1m1realp[i]+l1m1realp[i+1])/2 \
                                for i in range(0,l1r-1)]
                        points.append(l1m1realp[l1r-1]*3/2)
                        signes = [sign(abs(m2(x = s).real())-1) for s in points]
                        lt1 = len([x for x in signes if x >=0]) == 0
            self.properties["is_module_of_stability_function_less_than_1"] = lt1
            return lt1
    def is_A_stable(self):
        """
        Do we have an A_stable method ?

        See: HW II , second edition, page 43.
        """
        if "is_A_stable" in self.properties:
            return self.properties["is_A_stable"]
        else:
            ret = (self.is_module_of_stability_function_constant_on_Im() or \
                    self.is_module_of_stability_function_less_than_1()) and  \
                    self.real_part_of_poles_all_positive()[0]
            self.properties["is_A_stable"] = ret
            return ret
    def is_stiffly_accurate(self):
        """
        Do we have a stiffly accurate method ?
        """
        if "is_stiffly_accurate" in self.properties:
            return self.properties["is_stiffly_accurate"]
        else:
            # Proposition 3.8 H-W. TII, page 45.
            A = self.A
            B = self.B
            s1 = self.s-1
            Ls =  self.is_A_stable() and \
                (len([1 for j in range(0,self.s) if A[s1,j] != B[j]]) == 0 or \
                 len([1 for i in range(0,self.s) if A[i,0] != B[0]])  == 0)
            self.properties["is_stiffly_accurate"] = Ls
            return Ls
    def is_L_stable(self):
        """
        Do we have a L_stable method ?
        """
        if "is_L_stable" in self.properties:
            return self.properties["is_L_stable"]
        else:
            if not self.is_A_stable():
                ok= False
            else:
                R=self.stability_function()
                if R.denominator().degree() > R.numerator().degree():
                    ok = True
                elif  R.denominator().degree() <= R.numerator().degree():
                    ok = False
            self.properties["is_L_stable"] = ok
            return ok
    def is_algebraically_stable(self):
        """
        Is the method algebraically stable ?
        """
        if "is_algebraically_stable" in self.properties:
            return self.properties["is_algebraically_stable"]
        else:
            if self.is_explicit() or len([s for s in self.B if s<0])!=0:
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
                          print i,j,B[i]*A[i,j]+B[j]*A[j,i]-B[i]*B[j]
                        if not As:
                            break
            
            self.properties["is_algebraically_stable"] = As
            return As
    def conserve_quadratic_invariants(self):
        """
        Documentation is in the name of this method!
        """
        if "conserve_quadratic_invariants" in  self.properties:
            return self.properties["conserve_quadratic_invariants"]
        else:
            if self.is_explicit():
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
            self.properties["conserve_quadratic_invariants"] = ret
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
    def order_using_rooted_trees(self):
        """
        Compute the order of the method using rooted trees.
        """
        o = 1
        while self.check_order_using_rooted_trees(o):
            o+= 1
        return o-1
    
    def star_function(self,x,y):
        """
        Compute the star function. This is for drawing the "star" associated
        to the formula.
        """
        Rs = self.stability_function()
        s = Rs(x+QQbar(I)*y)/exp(x+I*y)
        return s*conjugate(s)
    def forget_properties(self):
        """
        Forget all what we already computed.
        """
        self.properties = {}
    def known_properties(self):
        """
        Return all the known properties at present time.
        """
        return self.properties
    def property(self,p):
        """
        Return the property 'p'
        """
        if  p in self.properties:
            return self.properties[p]
        else:
            return "None"
    def print_known_properties(self):
        for key in self.properties:
            print key,":", self.properties[key],"\n"
        
