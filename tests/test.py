#
import pytest
from rkkit import *
from methods.formulas import *
from rkkit.RKcolloc import *
#
class TestRKkit:
    """ Test computation of all known properties for a set of known
    classical methods """
    def testSetOfRKformulas(self):
        """ All known properties mus be computed with all methods.formulas."""
        cformulas=[Radau5(),Gauss4(),Radau2a(),RK4(),
                   SDIRK3(),SDIRK5(),Lobatto4()]
        Bads=[]
        #
        for f in cformulas:
    
            try:
                F=RKformula(f)
                F.compute_all_properties()
            except:
                Bads.append(f.Title)
        assert len(Bads) == 0
