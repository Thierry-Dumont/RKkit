#
import pytest
import passagemath_flint
from sage.rings.qqbar import *
import sage.rings.qqbar
#
from RKkit import *

#
class TestRKkit:
    """ Test computation of all known properties for a set of known
    classical methods """
    def testSetOfRKformulas(self):
        """ All known properties must be computed with all methods."""
        cformulas=[formulas.Radau5(),formulas.Gauss4(),formulas.Radau2a(),
                   formulas.RK4(),
                   formulas.SDIRK3(),formulas.SDIRK5(),formulas.Lobatto4()]
        Bads=[]
        #
        for f in cformulas:
            try:
                F=RKformula(f)
                print(f.Title)
                F.compute_all_properties()
            except:
                Bads.append(f.Title)
        assert len(Bads) == 0

