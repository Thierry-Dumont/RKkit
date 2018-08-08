

# Playing with Runge-Kutta methods and SageMath. #

_Beware: this is a very preliminary version!_

To use the codes, you must install [SageMath](http://www.sagemath.org/).

The goal of these programs is to verify different properties of a given Runge--Kutta method, defined by it's Butcher array (all the notions used here are defined in the _Bible_, see below).

* With SageMath, it is possible  to make *exact* computations in somme usefull sets of numbers, and thus we make the hypothesis that the coefficients of the methods live in the set of (real) algebraic numbers (called *AA* in SageMath), which is not a really restrictive hypothesis. Thus, the result obtained are  *proofs* (if my codes and SageMath are both correct!).

* To compute the order of a Runge-Kutta method, one must use the so called _rooted_ _trees_ for which we have a SageMath implementation, coded by [Florent Hivert](http://doc.sagemath.org/html/en/reference/combinat/sage/combinat/rooted_tree.html).

We also provide a function to compute the Butcher array of a method defined by collocation. A classical application is the set of Gaussian Runge-Kutta methods.

I also hope to implement the B-series in the future.

Two sage jupyter notebooks are provided: they probably provide  the best way to test this code.

#### References: ####

To learn about Ordinary Differential Equations solvers, you should read the
_Bible_:

*    Solving Ordinary Differential Equations I, by Hairer, Nørsett,, Wanner,
*    Solving Ordinary Differential Equations II Stiff and Differential-Algebraic
    Problems by Hairer and Wanner,
*    Geometric Numerical Integration by Hairer, Lubich and Wanner.


##### Note: ####
If you want to learn _SageMath_, you can read the book _Mathematical Computation
with Sage_ (which now is available in French, English and German), and
for which freely available pdf files are [downloadable 
here](https://members.loria.fr/PZimmermann/sagebook/english.html) and [there](http://sagebook.gforge.inria.fr/).


 # Using the code #



* First, you must define a Runge--Kutta method. To do this you must write a (simple) python class: have a look at one of "*.sage" files. 

This class must define:

  1- The arrays A and B of the Butcher array (the C part is not necessary).
  
  2- A title
  
  3- A comment.

Be carefull: The coefficients of A and B *must* live in algebraic numbers (AA or QQbar)!

* Then the best is to look at the notebook _Example1.ipynb_. For this, launch sage like this:

`>sage -n jupyter`

and launch the notebook _Exemple1.ipynb_.

#### Some cautions: ####

*  The code uses Sage's [@lazy_attribute](http://doc.sagemath.org/html/en/reference/misc/sage/misc/lazy_attribute.html) decorator; this allow to compute each
property only _once_ (as it can be very expensive).

This can be a bit disturbing:
 you must *not* do:

`sage: F.is_A_stable()` 

but:

`sage: F.is_A_stable`

as, with the @lazy_attribute decorator, refering to an attribute triggers the corresponding method ([see here](http://doc.sagemath.org/html/en/reference/misc/sage/misc/lazy_attribute.html)).

_Hint:_ if it does not work with brackets, try without, and vice-versa :-) .


## Gaussian formulae ##

Gausian formulae with n steps are obtained by collocation at the roots of the Legendre P polynomials of degree n, shifted from [-1,1] to [0,1].

RKcolloc.colloc computes the Butcher arrays (A and B parts);

See the notebook _Gaussian.ipynb._



### Notebooks ###

All are Sage / Jupyter notebooks. Launch Sage by typing:

`sage -n jupyter`


1. _Example.ipynb_ :  a tour of the system.

2. _Gaussian.ipynb_:  construct and test the Gauss methods.

3. _All properties.ipynb_: show how to compute all possible properties of a method.
