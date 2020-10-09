

# Playing with Runge-Kutta methods and SageMath. #

To use the code, you must install [SageMath](http://www.sagemath.org/).

The goal of this code is to verify different properties of a given Runge--Kutta method, defined by it's Butcher array (all the notions used here are defined in the _Bible_, see below).

* With SageMath, it is possible  to make *exact* computations in somme usefull sets of numbers, and thus we make the hypothesis that the coefficients of the methods live in the set of (real) algebraic numbers (called *AA* in SageMath); this is not a really restrictive hypothesis. But then, the result obtained are  *proofs* (if my codes and SageMath are both correct!).

* To compute the order of a Runge-Kutta method, one use the so called _rooted_ _trees_ for which we have a SageMath implementation, coded by [Florent Hivert](http://doc.sagemath.org/html/en/reference/combinat/sage/combinat/rooted_tree.html).

* We also provide a function to compute the Butcher array of a method defined by collocation. A classical application is the set of Gaussian Runge-Kutta methods.

I also hope to implement  B-series in the future.

Some sage/jupyter notebooks are provided: they probably provide  the best way to learn and to test this code.

#### References: ####

* To learn about Ordinary Differential Equations solvers, you should read the
_Bible_:

1.   Solving Ordinary Differential Equations I, by Hairer, Nørsett,, Wanner (HNW),
2.   Solving Ordinary Differential Equations II Stiff and Differential-Algebraic
         Problems by Hairer and Wanner 5HW),
3.   Geometric Numerical Integration by Hairer, Lubich and Wanner (HLW).

_If you want to learn and to understand rooted trees and B-series, I recommend to start first by reading reference 3.. Being the latest book, things have become much easier to understand than in 1._  



* If you want to learn _SageMath_, you can read the book _Mathematical Computation
with Sage_ (which now is available in French, English and German), and
for which freely available pdf files can be [downloaded here](https://members.loria.fr/PZimmermann/sagebook/english.html) and [there](http://sagebook.gforge.inria.fr/).


 # Using the code #

* First, you must define a Runge--Kutta method. To do this you must write a (simple) python _class_: have a look at one of the "*.sage" files in methods/.

This class must derive from "RungeKutta" and  include a constructor.

The constructor must:

- A) define:

  1- The arrays A and B of the Butcher array (the C part is not necessary).
  
  2- A title.

- B)  call the  base RungeKutta class contructor.

The "*.sage" files in methods/ give some examples of such classes.

Be carefull: The coefficients of A and B *must* be algebraic numbers (AA or QQbar) or rational numbers (QQ).

* Then the best is to look at the notebook _Example1.ipynb_. For this, launch sage like this:

`>sage -n jupyter`

and launch the notebook _Exemple1.ipynb_.

#### Implementation: ####

* The code uses a decorator @_persistance to avoid recomputing known properties (which is often expensive).

A precedent version was using  Sage's
[@lazy_attribute](http://doc.sagemath.org/html/en/reference/misc/sage/misc/lazy_attribute.html) decorator, which could be disturbing.


## Gaussian formulas (and other methods obtained be collocation) ##

Gausian formulae with n steps are obtained by collocation at the roots of the Legendre P polynomials of degree n, shifted from [-1,1] to [0,1].

RKcolloc.colloc computes the Butcher arrays, given collocation points in [0,1] and returns a Runge-Kutta method class  (note that collocation points are not necessary Gaussian points).

See the notebook _Gaussian.ipynb._



### Notebooks ###

All are Sage/Jupyter notebooks. Launch Sage by typing:

`sage -n jupyter`


1. _Example.ipynb_ :  a tour of the system.

2. _Gaussian.ipynb_:  construct and test some of the Gauss methods.

3. _AllProperties.ipynb_: show how to compute all possible properties of a method.

4. _Test.ipynb_ : just testing that everything works fine on a set of formulas.

5. _RadauByCollocation.ipynb_ : Radau methods are collocation methods. Play with this.
