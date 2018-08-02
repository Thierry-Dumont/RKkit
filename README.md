
### This repository will (soon) be populated with SageMath tools to play with Runge--Kutta  methods. ###

To use the codes, you must install [SageMath](http://www.sagemath.org/).

The goal of these programs is too verify the different properties of a given Runge--Kutta method, defined by it's Butcher array (all the names used here are defined in the _bible_, see below).

We use the possibily given by SageMath to make exact computations in interesting sets of numbers, and we make the hypothesis that the coefficients of the methods live in the set of (real) algebraic numbers (called *AA* in SageMath). Thus, the result obtained are  *proofs* (if my codes are correct!).

To compute the order of a Runge-Kutta method, one must use the so called _rooted_ _trees_ and SageMath provides an implementation of them, done by Florent Hivert, which we use here.

We also provide a function to compute the Butcher array of a method defined by collocation points.

I also hope to implement the B-series in the future.


#### References: ####

To learn about Ordinary Differential Equations solvers, you should read the
_bible_:

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



* First, you must define a Runge--Kutta method. Have a look at one of ".sage"* files. This class must define:

  1- The arrays A and B of the Butcher array (the C part is not necessary).
  
  2- A title
  
  3- A comment.

Be carefull: The coefficients of A and B *must* live in algebraic numbers (AA or QQbar)!

* Then the best is to look at the notebook _Example1.ipynb_. For this, launch sage like this

`> sage -n jupyter`

and lauch the notebook.
