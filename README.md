Stepalize
=========

Stepalize is a small, standalone Matlab program for the identification of
linear, time-invariant (LTI) dynamic models from step responses. It has been
successfully applied to the identification thermal dynamics for electronics
packaging, though it has many potential uses.


Features
--------
- Generates a linear, discrete-time estimate of a dynamic system from a
  measured step response.
    - Works with multi-output step responses.
    - Allows for model order selection via GUI at runtime.
- With YALMIP and SDPT3 installed, the following convex constraints are also
  possible:
    - Fix eigenvalues to lie within convex regions of the complex plane
          (such as the unit circle, real number line, etc.).
    - Fix the steady-state value of the system estimate.
    - Constrain estimates to have no overshoot and/or no undershoot.

The method uses a variant of the Ho-Kalman-Kung algorithm for impulse
responses but generalized to step responses. This is more accurate and far
more robust than trying to estimate an impulse response from differentiating
a step response.

In addition to identifying LTI systems, the method can be used to optionally
place constraints on the eigenvalues of the identified model. The
eigenvalues may be constrained to lie within any convex region of the
complex plane that is an intersection of ellipses, parabolas, and half
spaces that is symmetric about the real axis. Constraints that prevent the
resulting model from having any overshoot or undershoot in its step response
are also feasible.

Applying the constraints requires the (free) semidefinite program solver
[`SDPT3`](http://www.math.nus.edu.sg/~mattohkc/sdpt3.html) and the modeling
language parser [`YALMIP`](http://users.isy.liu.se/johanl/yalmip/).


Installation
------------
Just make sure `stepalize.m` is in your path. If you intend to use
semidefinite constraints, then make sure `SDPT3`, and `YALMIP` are in your
path and run `stepalize_test.m` to make sure everything works.

**NOTE**: The semidefinite program solver sometimes eats up a lot of memory,
and performance will vary quite a bit based on computing power. Try a small
problem first to make sure nothing crashes.


Documentation
-------------
See the help documentation for usage. Theoretical discussion can be found in
a preprint of the paper "Thermal Dynamical Identification of LEDs by
Step-Based Realization and Convex Optimization," which will appear in *IEEE
Transactions on Components, Packaging, and Manufacturing Technology*
sometime in 2013. It will have the doi
[10.1109/TCPMT.2012.2229464](http://dx.doi.org/10.1109/TCPMT.2012.2229464)
once it is published. (The link should work in the future.) A preprint may
be found
[here](https://sites.google.com/site/dnmiller/2012_Miller_IEEE_TCPMT.pdf).

