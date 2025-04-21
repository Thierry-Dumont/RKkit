from sage.all import *
from .RKExceptions import GraphicProblem
from sage.rings.infinity import minus_infinity

def RKplot(RKf,title="",Enlarge=4,TranslateX=0,
           ncurves=1,limits=[],fill=False,type="stab"):
    r"""
    Plot isovalues of stability function or Order star.

    PARAMETERS:
    ----------

    RKf: the formula (instantiation of RKformula).

    Enlarge: plot is done in a window around (0,0) in \mathbb{C}.

             We can Enlarge the size of the window by this factor (a real
             positive number).

    TranslateX : translate the origin in the window along x axis.

    ncurves: number of curves to plot.

    limits: [[min of x, max of x],[min of y, max of y]] for the window.
            In most case, this is computed.

    fill: to fill the plot (see contour_plot documentation).

    type= 'stab' for stability function (default), 'star' for the Order star.

    """
    def sf1(x,y,P):
        s = P(x+QQbar(I)*y)
        return (s*conjugate(s)).real()

    RDroots = RKf.poles_of_stability_function()
    Rstab= RKf.stability_on_real_negative_axis()

    # try to compute limits:
    if Rstab==minus_infinity and RDroots == [] :
        if limits==[]:
            raise GraphicProblem(
                "Cannot find limits. You must define: limits = [(x1,y1),(x2,y2)]")
        else:
            raise GraphicProblem("Mthod is not A_stable and no stab. limit is known")
          
    else:
        if limits == []:
            if len(RDroots) > 0:
                Lr = max(max([abs(s[0].n().real()) for s in RDroots]), \
                        max([abs(s[0].n().imag()) for s in RDroots]))
            else:
                Lr = 0
                
            if Rstab != minus_infinity:
                L = max(-Rstab,Lr)
                Lm= -L*Enlarge
                Lp= -Lm
            elif RKf.is_A_stable:
                L = Lr
                Lm = -Enlarge*L
                Lp = 3*Enlarge*L/2
            else:
                Lm = -L
                Lp = L
            
            L1 = (Lm,Lp)
            L2 = (-Enlarge*L,Enlarge*L)
            
        #
    if limits != []:
        L1 = limits[0]
        L2 = limits[1]

    # What to plot:
    if type ==   "stab":
        sf = lambda x,y: sf1(x,y,RKf.stability_function())
    elif type == "star":
        sf = RKf.order_star_function()
    else:
        raise AttributeError("RKplot: "+type+" :unknown plot type")

    if title != "":
        stitle =  title
    else:
        if type ==  "stab":
            stitle= "Stability function."
        else:
            stitle = "Order star."
            
    # translate:
    if TranslateX!=0:
        d = (L1[1]-L1[0])*TranslateX/100.
        L1 = (L1[0]-d,L1[1]-d)
        
    x = SR.var("x")
    y = SR.var("y")

    return contour_plot(sf,(x,L1[0],L1[1]), (y,L2[0],L2[1]),
                     contours = [0,1],
                     labels = True,fill = fill, label_inline = True,
                     axes = True,colorbar = True,title = stitle)                 

