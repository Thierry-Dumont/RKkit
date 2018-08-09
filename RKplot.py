from __future__ import print_function
from sage.symbolic.ring import SR
from sage.symbolic.constants import I
from sage.rings.all import QQbar
from sage.functions.other import conjugate
from sage.plot.contour_plot import *
from sage.plot.point import *
from RKExceptions import GraphicProblem
from sage.rings.infinity import minus_infinity

def RKplot(RKf,title="",Enlarge=4,ncurves=1,limits=[],fill=False,type="stab"):
    def sf1(x,y,P):
        s = P(x+QQbar(I)*y)
        return s*conjugate(s)

    RDroots = RKf.poles_of_stability_function
    Rstab= RKf.stability_on_real_negative_axis

    if Rstab==minus_infinity and RDroots == [] :
        if limits==[]:
            raise GraphicProblem(
                "Cannot find limits. You must define: limits = [(x1,y1),(x2,y2)]")
        else:
            raise GraphicProblem("Mthod is not A_stable and no stab. limit is known")
            #L1 = limits[0]
            #L2 = limits[1]
    else:
        if limits == []:
            if len(RDroots)>0:
                Lr= max(max([abs(s[0].n().real()) for s in RDroots]), \
                        max([abs(s[0].n().imag()) for s in RDroots]))
            else:
                Lr=0
                
            if Rstab != minus_infinity:
                L=max(-Rstab,Lr)
                Lm= -L*Enlarge
                Lp= -Lm
            elif RKf.is_A_stable:
                L= Lr
                Lm = -Enlarge*L
                Lp = 3*Enlarge*L/2
            else:
                Lm=-L
                Lp = L
            
            L1 = (Lm,Lp)
            L2 = (-Enlarge*L,Enlarge*L)
        #
    if limits != []:
        L1 = limits[0]
        L2 = limits[1]
    if type ==   "stab":
        sf=lambda x,y: sf1(x,y,RKf.stability_function)
    elif type == "star":
         sf=lambda x,y: RKf.order_star_function(x=x,y=y)
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
    x=SR.var("x")
    y=SR.var("y")
    if withCmap:
        c = contour_plot(sf,(x,L1[0],L1[1]), (y,L2[0],L2[1]),
                         contours=[-1,-0.5,0,1],
                         labels=True,fill=fill, label_inline=True,axes=True,colorbar=True)
                       
    else:
        c = contour_plot(sf,(x,L1[0],L1[1]), (y,L2[0],L2[1]), \
                         contours=[-1,0,1],labels=True,label_inline=True,
                         axes=True,fill=fill,colorbar=True)
                       
    return c

