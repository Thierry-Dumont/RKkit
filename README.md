

# Playing with Runge-Kutta methods and SageMath. #

_Beware: this is a very preliminary version!_

To use the code, you must install [SageMath](http://www.sagemath.org/).

The goal of this code is to verify different properties of a given Runge--Kutta method, defined by it's Butcher array (all the notions used here are defined in the _Bible_, see below).

* With SageMath, it is possible  to make *exact* computations in somme usefull sets of numbers, and thus we make the hypothesis that the coefficients of the methods live in the set of (real) algebraic numbers (called *AA* in SageMath); this is not a really restrictive hypothesis. But then, the result obtained are  *proofs* (if my codes and SageMath are both correct!).

* To compute the order of a Runge-Kutta method, one must use the so called _rooted_ _trees_ for which we have a SageMath implementation, coded by [Florent Hivert](http://doc.sagemath.org/html/en/reference/combinat/sage/combinat/rooted_tree.html).

* We also provide a function to compute the Butcher array of a method defined by collocation. A classical application is the set of Gaussian Runge-Kutta methods.

I also hope to implement  B-series in the future.

Some sage/jupyter notebooks are provided: they probably provide  the best way to learn and to test this code.

#### References: ####

* To learn about Ordinary Differential Equations solvers, you should read the
_Bible_:

1.    Solving Ordinary Differential Equations I, by Hairer, Nørsett,, Wanner (HNW),
2.   Solving Ordinary Differential Equations II Stiff and Differential-Algebraic
    Problems by Hairer and Wanner 5HW),
3.   Geometric Numerical Integration by Hairer, Lubich and Wanner (HLW).

_If you want to learn and to understand rooted trees and B-series, I recommend to start first by reading reference 3.. Being the latest book, things have become much easier to understand than in 2._  



* If you want to learn _SageMath_, you can read the book _Mathematical Computation
with Sage_ (which now is available in French, English and German), and
for which freely available pdf files can be [downloaded](https://members.loria.fr/PZimmermann/sagebook/english.html) and [there](http://sagebook.gforge.inria.fr/).


 # Using the code #

* First, you must define a Runge--Kutta method. To do this you must write a (simple) python class: have a look at one of "*.sage" files. 

This class must define:

  1- The arrays A and B of the Butcher array (the C part is not necessary).
  
  2- A title
  

Be carefull: The coefficients of A and B *must* be algebraic numbers (AA or QQbar) or rational numbers (QQ).

* Then the best is to look at the notebook _Example1.ipynb_. For this, launch sage like this:

`>sage -n jupyter`

and launch the notebook _Exemple1.ipynb_.

#### Some cautions: ####

*  The code uses Sage's [@lazy_attribute](http://doc.sagemath.org/html/en/reference/misc/sage/misc/lazy_attribute.html) decorator; this allow to compute _each_
_property_ _only_ _once_ (as computations  can be very expensive).

This can be a bit disturbing:
 do *not* write:

`sage: F.is_A_stable()` 

but write:

`sage: F.is_A_stable`

as, with the @lazy_attribute decorator, refering to an attribute triggers (only once!) the corresponding method ([see here](http://doc.sagemath.org/html/en/reference/misc/sage/misc/lazy_attribute.html)).

Only two methods must be called as usual, as they do not create attributes:

`sage: F.compute_all_properties()`

and:

`sage: F.print_all_known_properties()`

_Hint:_ if something does not work with brackets, try without, and vice-versa :smile:.


## Gaussian formulae ##

Gausian formulae with n steps are obtained by collocation at the roots of the Legendre P polynomials of degree n, shifted from [-1,1] to [0,1].

RKcolloc.colloc computes the Butcher arrays (A and B parts);

See the notebook _Gaussian.ipynb._



### Notebooks ###

All are Sage/Jupyter notebooks. Launch Sage by typing:

`sage -n jupyter`


1. _Example.ipynb_ :  a tour of the system.

2. _Gaussian.ipynb_:  construct and test the Gauss methods.

3. _All properties.ipynb_: show how to compute all possible properties of a method.

4. _Test.ipynb_ : just testing that all work fine on a set of formulas.
