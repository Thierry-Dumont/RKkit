# -*- coding: utf-8 -*-
from  sage.all import *
from sage.combinat.rooted_tree import RootedTree as RT
from sage.combinat.rooted_tree import RootedTrees_size as RTS
#
class  RKTrees(SageObject):
    r"""
    The rooted trees machinery.
    
    EXAMPLES:

    sage: R= RKTrees(n)

    n is the maximum depth of the rooted trees you will use (the
    built dictionary will be enlarged if necessary (lazzy evaluation)).

    Bibliography: HW are the books of Hairer, Wanner and co-workers.
    """
    def __init__(self):
        self.n = 1
        self.dtrees = {}
        self.expand(1)
    def expand(self,l):
        r"""
        Extend the list of rooted trees up to depth l
        """
        for i in range(1,l+1):
            if not i in self.dtrees:
                self.dtrees[i] = RTS(i).list()
        self.n=l
    def gamma(self,t):
        r"""
        The gamma coefficient (see HW referenced books).
        """
        tn = t.node_number()
        if tn == 1:
            return 1
        else:
            return tn*prod([self.gamma(s) for s in t])
    def _LabelledTree_to_formula(self,rtc,root_label,faclist):
        if rtc.node_number() == 1:
            faclist.append((root_label,rtc.label()))
        else:
            for t in rtc:
                faclist.append((root_label,t.label()))
                if t.node_number() > 1:
                    self._LabelledTree_to_formula(t,t.label(),faclist)
    def tree_to_order_formula(self,rt):
        rtc = rt.canonical_labelling()
        faclist = []
        self._LabelledTree_to_formula(rtc,rtc.label(),faclist)
        return faclist
    def eval_sum_prod(self,A,B,formula,v):
        return  B[v[0]] * prod( [A[v[i[0]-1],v[i[1]-1]] for i in formula] )
    def tree_order_form(self,A,B,rt):
        n = len(B)
        s = set(range(0,n))
        S = cartesian_product([s for i in range(0,rt.node_number())])
        #
        f = self.tree_to_order_formula(rt)
        #
        s =  sum([self.eval_sum_prod(A,B,f,v) for v in S])
        return s*self.gamma(rt) 
    def check_tree_order(self,A,B,rt):
        return self.tree_order_form(A,B,rt) == 1
    def check_order(self,A,B,order):
        """
        Prove that the rooted trees of order 'order' fullfill 
        the requirements.

        To check if a Runge-Kutta formula as order 'n', one must call 
        check_order(A,B,i) for i in range(1,n+1).
        """
        if order == 1:
            s = sum(B)
            s.exactify()#should make sometime more readable results, may be.
            return s == 1
        else:
            for i in range(len(self.dtrees)+1,order+1):
                self.expand(i)
            for t in self.dtrees[order]:
                tc = t.canonical_labelling()
                ok = self.check_tree_order(A,B,tc)
                if not ok:
                    return False
            return True
 
    def symetry_coefficient(self,rt):
        rt1 = RT(rt)
        if rt1 ==  RT([]):
            return 1
        else:
            l = [t for t in rt1]
            l.sort()
            ft = 1
            f = 1
            for i in range(0,len(l)-1):
                if l[i] == l[i+1]:
                    f+= 1
                else:
                    ft*= factorial(f)
                    f=1
            if f != 1: ft*= factorial(f)
            return prod([self.symetry_coefficient(s) for s in l])*ft
    def compute_gamma(self):
        r"""
        This returns the '\gamma' coefficients (see HW).
        """
        self.gamma={}
        for i in range(1,self.n+1):
            for s in self.dtrees[i]:
                self.gamma[s] = self.symetry_coefficient(s)
