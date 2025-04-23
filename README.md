

# Playing with Runge-Kutta methods and SageMath. #

To use the code, you must install [SageMath](http://www.sagemath.org/).

The goal of this code is to verify different properties of a given Runge--Kutta method, defined by it's Butcher array (all the notions used here are defined in the _Bible_, see below).

* With SageMath, it is possible  to make *exact* computations in some useful sets of numbers, and thus we make the hypothesis that the coefficients of the methods live in the set of (real) algebraic numbers (called *AA* in SageMath); this is not a really restrictive hypothesis. But then, the result obtained are  *proofs* (if my codes and SageMath are both correct!).

* To compute the order of a Runge-Kutta method, one use the so called _rooted_ _trees_ for which we have a SageMath implementation, coded by [Florent Hivert](http://doc.sagemath.org/html/en/reference/combinat/sage/combinat/rooted_tree.html).

* We also provide a function to compute the Butcher array of a method defined by collocation. A classical application is the set of Gaussian Runge-Kutta methods.

Some sage/jupyter notebooks are provided: they probably provide  the best way to learn and to test this code (even, you can run the notebooks on binder, see below).



#### References: ####

* To learn about Ordinary Differential Equations solvers, you should read the
_Bible_:

1-   Solving Ordinary Differential Equations I, by Hairer, Nørsett,,
Wanner (HNW),

2-   Solving Ordinary Differential Equations II Stiff and Differential-Algebraic
         Problems by Hairer and Wanner (HW),
		 
3-   Geometric Numerical Integration by Hairer, Lubich and Wanner (HLW).

_If you want to learn and to understand rooted trees and B-series, I
recommend to start first by reading reference 3. Being the latest
book, things have become much easier to understand than in reference
1._   

* If you want to learn _SageMath_, you can read the book _Mathematical Computation
with Sage_ (which now is available in French, English and German), and
for which freely available pdf files can be [downloaded
here](https://members.loria.fr/PZimmermann/sagebook/english.html) and
[there](http://sagebook.gforge.inria.fr/). 


## _Using the code_ 

* __First__, you must define a Runge--Kutta method. To do this you
  must write a (simple) python _class_: have a look at _formulas.py_
  in methods/ to undestand what to do.

	The class must derive from the "RungeKutta" class and  include a
    constructor. 

The constructor must:


1- Define the arrays A and B of the Butcher array (the C part is generally
    not necessary). 
  
2- - Give a title.

3-  Call the  base RungeKutta class constructor.

Remember that the  "formula.py" file in methods/ give  examples of such classes.

__Be careful:__ The coefficients of A and B *must* be algebraic numbers
(AA or QQbar) or rational numbers (QQ). This is not absolutely trivial
since, in sage:

2/3, for example, is directly evaluated as a float (0.66666...). This
is a Python intrinsic, and cannot be avoided. So as you need to enter
an algebraic (real) number or a rational, you _must_ write QQ(2/3) or
AA(2/3) or 2/QQ(3) or 2/AA(3).

Have a look at "GoodAndBad.ipynb": this notebook shows this more in details.


* __Then__, the best is to look at the notebook _Example1.ipynb_ and
  adapt it to your needs. For this, launch sage like this:

`>sage -n jupyter`

or

`>sage -n jupyterlab` 

and then, launch the notebook _Exemple1.ipynb_.



## Gaussian formulas (and other methods obtained by collocation) ##

Gausian formulae with n steps are obtained by collocation at the roots of the Legendre P polynomials of degree n, shifted from [-1,1] to [0,1].

RKcolloc.colloc computes the Butcher arrays, given collocation points in [0,1] and returns a Runge-Kutta method class  (note that collocation points are not necessary Gaussian points).

See the notebook _Gaussian.ipynb._



### Notebooks ###

All are Sage/Jupyter notebooks. Launch Sage by typing:

`sage -n jupyter`

or

`sage -n jupyterlab`

1. _Example.ipynb_ :  a tour of the system.

2. _Gaussian.ipynb_:  construct and test some of the Gauss methods.

3. _AllProperties.ipynb_: show how to compute all possible properties of a method.

4. _Test.ipynb_ : just testing that everything works fine on a set of formulas.

5. _RadauByCollocation.ipynb_ : Radau methods are collocation methods. Play with this.

6. _GoodAndBad.ipynb_ : is supposed to show what to do and what not to
   do when coding a Runge-Kutta formula.

### Running the notebooks on Binder ###
Just 
[Click me.](https://mybinder.org/v2/gh/Thierry-Dumont/RKkit/315376e77071abff5ab16ab9f6ecba52a3c359e0)

otherwise, if it does not work, have a look at  "Dockerfile".

#### Implementation: ####

* The code uses a decorator @_persistance to avoid recomputing known properties (which is often expensive).

A precedent version was using  Sage's
[@lazy_attribute](http://doc.sagemath.org/html/en/reference/misc/sage/misc/lazy_attribute.html) decorator, which could be disturbing.
