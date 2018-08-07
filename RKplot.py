from sage.symbolic.ring import SR
from sage.symbolic.constants import I
from sage.rings.all import QQbar
from sage.functions.other import conjugate
from sage.plot.contour_plot import *
from sage.plot.point import *
from RKExceptions import GraphicProblem

def RKplot(RKf,title="",Enlarge=4,ncurves=1,limits=[],fill=False,type="stab"):
    def sf1(x,y,P):
        s = P(x+QQbar(I)*y)
        return s*conjugate(s)
    def sfStar(x,y,F):
        return F.star_function(x,y)
    #sfStar = lambda x,y:  F.star_function(x,y)
    #sf1 =  lambda x,y:  P(x+QQbar(I)*y)
    RDroots = RKf.poles_of_stability_function()
    if RKf.is_explicit() or  RDroots == [] :
        if limits==[]:
            raise GraphicProblem(
                "Cannot find limits. You must define: limits = [(x1,y1),(x2,y2)]")
        else:
            L1 = limits[0]
            L2 = limits[1]
    else:
        if limits == []:
            L = Enlarge*max(max([abs(s[0].n().real()) for s in RDroots]), \
                          max([abs(s[0].n().imag()) for s in RDroots]))

         
            if RKf.is_A_stable():
                Lm = 0
                Lp = 2*L
            else:
                Lm=-L
                Lp = L
            L1 = (Lm,Lp)
            L2 = (-L,L)
        #
    if limits != []:
        L1 = limits[0]
        L2 = limits[1]
    if type ==   "stab":
        sf=lambda x,y: sf1(x,y,RKf.stability_function())
    elif type == "star":
        sf =  lambda x,y: sfStar(x,y,RKf)
    else:
        raise AttributeError("RKplot: "+type+" :unknown plot type")
    if ncurves == 1:
        contrs = [1]
        withCmap =  True
    else:
        h = 1./(ncurves-1)
        contrs = [1+h*i for i in range(0,ncurves)]
        withCmap = False
    if fill: contrs.append(99999)

    if title != "":
        stitle =  "Stability domain for "+title
    else:
        stitle = ""
    if withCmap:
        c = contour_plot(sf,L1,L2,fill=fill, \
                       contours=contrs,cmap=[(0,1,0),(1,0,0)], \
                       labels=True,title=stitle)
    else:
        c = contour_plot(sf,L1,L2,fill=fill, \
                       contours=contrs, \
                       labels=True,title=stitle)
    return c

