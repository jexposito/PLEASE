/*
 * optimpack.i -
 *
 *	Optimization routines for Yorick (OptimPack version 2).
 *
 *
 * Acknowledgments:
 *	This package is partially based on routines written by others:
 *
 *	 - The line search and variable metric limited memory routines
 *	   are derived from the MINPACK-2 project:
 *           ftp://info.mcs.anl.gov/pub/MINPACK-2/
 *
 *	 - Conjugate gradient is based on CG+ package by Gilbert, J.C. and
 *	   Nocedal, J. [1]: http://www.ece.nwu.edu/~rwaltz/CG+.html
 *
 *	The original routines were converted in Yorick and _improved_ in
 *	many aspects:
 *
 *	 1. The multidimensional minimization routines automatically
 *	    resume the search whenever the computed search direction is
 *	    not a descent and should therefore be more robust with respect
 *	    to non-linear problems.
 *
 *	 2. Preconditioning can be used in conjugate gradient method.
 *
 *       3. The multidimensional minimization routines can account for
 *	    simple bound constraints on the parameters.
 *
 *
 * References:
 *	[1] Gilbert, J.C. and Nocedal, J. (1992). "Global Convergence
 *	    Properties of Conjugate Gradient Methods", SIAM Journal on
 *	    Optimization, Vol. 2, pp. 21-42.
 *	[2] Shewchuk, J.R. (1994). "An Introduction to the Conjugate
 *	    Gradient Method Without the Agonizing Pain",
 *	    ftp://warp.cs.cmu.edu/quake-papers/painless-conjugate-gradient.ps
 *      [3] A. Schwartz and E. Polak, "Family of Projected Descent Methods for
 *	    Optimization  Problems with Simple Bounds", Journal of
 *	    Optimization Theory and Applications, Vol. 92, pp. 1-31, 1997.
 *
 *-----------------------------------------------------------------------------
 *
 *	Copyright (c) 2001-2003 Eric THIEBAUT.
 *
 *	This file is part of OptimPack.
 *
 *	OptimPack is  free software; you can redistribute  it and/or modify
 *	it under the  terms of the GNU General  Public License as published
 *	by the Free  Software Foundation; either version 2  of the License,
 *	or (at your option) any later version.
 *
 *	OptimPack is  distributed in the hope  that it will  be useful, but
 *	WITHOUT  ANY  WARRANTY;  without   even  the  implied  warranty  of
 *	MERCHANTABILITY or  FITNESS FOR A PARTICULAR PURPOSE.   See the GNU
 *	General Public License for more details.
 *
 *	You should have  received a copy of the  GNU General Public License
 *	along with  Yeti (file "COPYING"  in the top source  directory); if
 *	not, write to the Free  Software Foundation, Inc., 59 Temple Place,
 *	Suite 330,  Boston, MA  02111-1307 USA.
 *
 *-----------------------------------------------------------------------------
 *
 * History:
 *	$Id: optimpack.i,v 1.1 2005/06/10 06:19:40 eric Exp eric $
 *	$Log: optimpack.i,v $
 *	Revision 1.1  2005/06/10 06:19:40  eric
 *	Initial revision
 *
 */

/*
 * CHANGES:
 *  - Always use effective step (option PRJ removed).
 *  - Routines renamed after those in C-version of OptimPack-2.
 *  - Use safeguarded line search step based on how far is the
 *    farest bound (option FMIN removed, new option SXBIG).
 *  - Use oldest (S,Y) pair in VMLM to store initial parameter and
 *    (projected) gradient to save memory.
 *
 * THINGS TO DO:
 *  - Fix Shanno & Phua in conjugate-gradient methods.
 *  - Projection in conjugate-gradient methods.
 *  - Fix computation of AMAXG at start.
 *  - Reduce number of computations of FREE and GP = FREE*G.
 *  - Fix documentation.
 *  - Fix the "rounding errors prevent progress" bug (may be a clever
 *    restart of the algorithm).
 *  - Use "best solution found so far".
 *  - Fix restart (should be done more often).
 */

/*---------------------------------------------------------------------------*/
/* Save built-in version, if any, of functions implemented by OptimPack
   (these statements _must_ be at the beginning of this file and are
   produced by the script optimpack-fix). */

if (is_func(op_cgmnb_next)==2 && is_func(op_cgmnb_next_BUILTIN)!=2)
  op_cgmnb_next_BUILTIN = op_cgmnb_next;
if (is_func(op_cgmnb_setup)==2 && is_func(op_cgmnb_setup_BUILTIN)!=2)
  op_cgmnb_setup_BUILTIN = op_cgmnb_setup;
if (is_func(op_get_amaxg)==2 && is_func(op_get_amaxg_BUILTIN)!=2)
  op_get_amaxg_BUILTIN = op_get_amaxg;
if (is_func(op_get_flags)==2 && is_func(op_get_flags_BUILTIN)!=2)
  op_get_flags_BUILTIN = op_get_flags;
if (is_func(op_get_iter)==2 && is_func(op_get_iter_BUILTIN)!=2)
  op_get_iter_BUILTIN = op_get_iter;
if (is_func(op_get_msg)==2 && is_func(op_get_msg_BUILTIN)!=2)
  op_get_msg_BUILTIN = op_get_msg;
if (is_func(op_get_rejects)==2 && is_func(op_get_rejects_BUILTIN)!=2)
  op_get_rejects_BUILTIN = op_get_rejects;
if (is_func(op_get_restarts)==2 && is_func(op_get_restarts_BUILTIN)!=2)
  op_get_restarts_BUILTIN = op_get_restarts;
if (is_func(op_get_stage)==2 && is_func(op_get_stage_BUILTIN)!=2)
  op_get_stage_BUILTIN = op_get_stage;
if (is_func(op_get_step)==2 && is_func(op_get_step_BUILTIN)!=2)
  op_get_step_BUILTIN = op_get_step;
if (is_func(op_vmlmb_next)==2 && is_func(op_vmlmb_next_BUILTIN)!=2)
  op_vmlmb_next_BUILTIN = op_vmlmb_next;
if (is_func(op_vmlmb_setup)==2 && is_func(op_vmlmb_setup_BUILTIN)!=2)
  op_vmlmb_setup_BUILTIN = op_vmlmb_setup;

/*---------------------------------------------------------------------------*/
local op_info;
/* DOCUMENT OptimPack - a Yorick package for optimization.

     Routines:
        op_fmin - minimization of an univariate function.
        op_mnb - minimization of an multivariate function, possibly
              with bound constraints.
        op_cgmnb_setup, op_cgmn_next - conjugate gradients algorithm.
        op_vmlmb_setup, op_vmlmb_next - VMLM-B algorithm.
 */

/*---------------------------------------------------------------------------*/
/* MINIMIZATION OF AN UNIVARIATE FUNCTION */

func op_fmin(f, a, b, lim, tol=, all=, eps=)
/* DOCUMENT op_fmin(f, a, b)
       -or- op_fmin(f, a, b, lim)
     Get the location  of the minimum of univariate function  F(X).  F is a
     Yorick  function,  A and  B  are bounds  or  starting  values for  the
     variable X and optional LIM specifies the kind of limits for X:
       If LIM=0 or nil, there is no bounds for X: F is first evaluated at A
         and B, then the interval of  search is enlarged until a minimum is
         bracketed.
       If LIM=1,  the interval is bounded by  A: F is first  evaluated at B
         and  the interval is  enlarged (away  from A)  until a  minimum is
         bracketed  --  i.e.  the location of  the  minimum X is such that:
         A < X, if A < B; or A > X, if A > B.
       If LIM=2, the  interval is bounded by B (same as  with LIM=1 but the
         role of A and B exchanged).
       If LIM=3, the minimum is searched  in the interval (A,B) -- i.e. the
         location of the minimum X is such that: min(A,B) < X < max(A,B).

     Keyword  TOL can be  used to  specify the  relative precision  for the
     solution.  The default value is TOL=sqrt(EPS) (see below).

     If keyword ALL  is true (non-nil and non-zero)  the returned value is:
     [X, FX, XLO, XHI] where X is the approximated location of the minimum,
     FX is F(X) and XLO and XHI are the lower and upper bounds for the true
     minimum; otherwise, only X is returned.

     Keyword EPS  can be  used to specify  the machine  relative precision.
     Default value is EPS~2.22e-16,  which corresponds to IEEE standard for
     double precision.

   NOTES:
     (1) The  minimization routine  should never evaluates  F at  the given
         bounds if any.
     (2) If the function F(X) is  not unimodal, only a local minimum can be
         found.

   REFERENCES:
     The method is  based on original Brent's method  modified to allow for
     different kind of bounds (both, left, right or none).

   SEE ALSO: op_mnb. */
{
  /* Make sure A and B are double precision values. */
  a += 0.0;
  b += 0.0;

  /* EPS is approximately the square root of the relative machine
     precision. */
  if (is_void(eps)) eps = 2.2204460492503131e-16; /* assume IEEE double */
  tol1 = eps + 1.0;
  eps = sqrt(eps);
  tol3 = (is_void(tol) ? eps : tol)/3.0;
  /* TOL not used below */

  /* C = (3 - sqrt(5))/2 is the squared inverse of the golden ratio */
  c = 0.3819660112501051517954131656343618822796908201942371378645513772947395;

  /* S = (1 + sqrt(5))/2 = 2 - C is a constant used to increase the width
     of the interval with a golden ratio. */
  s = 1.6180339887498948482045868343656381177203091798057628621354486227052605;

  /* Original Brent's method assumes that the minimum is in (A,B) with
   * A<=B and keeps track of the following variables:
   *   X, FX = least function value found so far
   *   W, FW = previous value of X, FX
   *   V, FV = previous value of W, FW
   *   U, FU = last function evaluation
   * If the interval to consider is not bounded or only left/right bounded,
   * the idea is to find a suitable interval (A,B) where at least one
   * minimum must exists (if the function is continue) and start Brent's
   * algorithm with correct values for X, FX, ... (in order to save some
   * function evaluations).
   */
  if (! lim) {
    /* The interval of search is unlimited, we start with A, B and then
       search for a bracket. */
    x = a;
    fx = f(x);
    w = b;
    fw = f(w);

    /* Make sure X is the best location found so far, we therefore exchange
       W and X if FW <= FX (we exchange the two point in case of equality
       to alternatively search on the other side where the function is,
       numerically, flat) */
    if (fw <= fx) {
      tmp=w; w=x; x=tmp;
      tmp=fw; fw=fx; fx=tmp;
    }

    /* Loop until a bracket is found. Possible improvements: (1) use parabolic
       extrapolation to allow for bigger jumps, (2) keep track of one more
       point in order to be able to slightly reduce the size of the interval
       once a bracket has been found. */
    for (;;) {
      /* Take a golden step in the descent direction. */
      v = x + s*(x - w);
      fv = f(v);

      if (fw > fx) {
        if (fv > fx) {
          /* Bracket found: the minimum is in (V,W).  Set variables for
             Brent's method: set bounds such that A is smaller than B. */
          if (v < w) {
            a = v;
            b = w;
          } else {
            a = w;
            b = v;
          }
          break; /* branch to Brent's method */
        } else {
          /* Continue with golden search with X and V. */
          w=x; fw=fx;
          x=v; fx=fv;
        }
      } else {
        /* We are in a, numerically, flat region (FW=FX). Enlarge interval
           (V,W) by a factor 1+S ~ 2.62 */
        if (fv >= fx) {
          x=w; fx=fw;
          w=v; fw=fv;
        } else {
          w=x; fw=fx;
          x=v; fx=fv;
        }
      }
    }
  } else if (lim == 1 || lim == 2) {
    /* Interval is bounded by A or B.  Possibly exchange A and B, so that
       A is the bound and search until a bound for the other side is found. */
    if (lim == 2) {tmp=a; a=b; b=tmp;}
    w = x = b;
    fw = fx = f(x);
    for (;;) {
      v = x + s*(x - a);
      fv = f(v);
      if (fv > fx) {
        /* We have found a bound for the other side.  Set search interval
           to be (A,B) := (A,V) or (V,A) such that A is smaller than B. */
        if (v > a) {
          b = v;
        } else {
          b = a;
          a = v;
        }
        break; /* branch to Brent's method */
      }
      w=x; fw=fx;
      x=v; fx=fv;
      if (fw > fx) a = w;
    }
  } else if (lim == 3) {
    /* The minimum is to be found in (A,B) -- this is original Brent's
       method.  Make sure that A is smaller than B and set X,FX ... */
    if (a > b) {tmp=a; a=b; b=tmp;}
    v = w = x = a + c*(b - a);
    fv = fw = fx = f(x);
  } else {
    error, "bad value for keyword LIM";
  }

  /*** Brent's method. ***/

  /* Set E and D (note: the golden step instead of the parabolic step is
     taken if abs(E) is too small). */ 
  e = x - v;
  d = x - w;

  /* main loop starts here */
  for (;;) {
    xm = (a + b)*0.5;
    tol1 = eps*abs(x) + tol3;
    tol2 = tol1 + tol1;

    /* check stopping criterion */
    if (abs(x - xm) <= tol2 - (b - a)*0.5) {
      if (all) return [x, fx, a, b];
      return x;
    }
    if (abs(e) > tol1) {
      /* fit parabola */
      q = (x - v)*(fx - fw);
      r = (x - w)*(fx - fv);
      if (q <= r) {
	p = (x - v)*q - (x - w)*r;
	q = (r - q)*2.0;
      } else {
	p = (x - w)*r - (x - v)*q;
	q = (q - r)*2.0;
      }
      if (abs(p) < abs(0.5*q*e) && p > q*(a - x) && p < q*(b - x)) {
	/* use a parabolic-interpolation step */
        e = d;
        u = x + (d = p/q);

	/* F must not be evaluated too close to A or B */
	if (u - a < tol2 || b - u < tol2) {
	  d = (x < xm ? tol1 : -tol1);
	}
      } else {
	/* use a golden-section step */
	e = (x >= xm ? a : b) - x;
	d = c*e;
      }
    } else {
      /* use a golden-section step */
      e = (x >= xm ? a : b) - x;
      d = c*e;
    }

    /* F must not be evaluated too close to X */
    u = (abs(d) >= tol1 ? x + d : (d > 0.0 ? x + tol1 : x - tol1));
    fu = f(u);

    /* update A, B, V, W, and X */
    if (fx <= fu) {
      if (u >= x) b = u;
      else        a = u;
    }
    if (fu <= fx) {
      if (u >= x) a = x;
      else        b = x;
      v = w; fv = fw;
      w = x; fw = fx;
      x = u; fx = fu;
    } else if (fu <= fw || w == x) {
      v = w; fv = fw;
      w = u; fw = fu;
    } else if (fu <= fv || v == x || v == w) {
      v = u; fv = fu;
    }
  } /* end of main loop */
}

/*---------------------------------------------------------------------------*/
/* SIMPLIFIED INTERFACE FOR LIMITED MEMORY MULTIDIMENSIONAL MINIMIZATION */

func op_z2d(z) { return [double(z), z.im]; }
func op_d2z(d) { return (1i*d(..,2) + d(..,1)); }
/* DOCUMENT d = op_z2d(z);
       -or- z = op_d2z(d);
     Convert between complex array Z and array of doubles D.  The real part
     of Z is D(..,1) and the imaginary part of Z is D(..,2). 

   SEE ALSO: double, complex. */

func op_mnb(f, x, &fout, &gout, &cpuout, &evalout, &steplenout, &maxgout,
            extra=, xmin=, xmax=, method=, mem=, maxiter=, maxeval=, verb=,
            ftol=, gtol=, sftol=, sgtol=, sxtol=, sxbig=)
/* DOCUMENT op_mnb(f, x)
       -or- op_mnb(f, x, fout, gout)
     Returns  a  minimum  of   a  multivariate  function  by  an  iterative
     minimization algorithm (conjugate  gradient or limited memory variable
     metric)  possibly with  simple  bound constraints  on the  parameters.
     Arguments are:

       F - User defined function to optimize.  The prototype of F is as
           follow:
             func F(x, &gx) {
               fx = ....; // compute function value at X
               gx = ....; // store gradient of F in GX
               return fx; // return F(X)
             }
       X - Starting solution (a floating point array).
       FOUT - Optional output variable to store the value of F at the
           minimum.           
       GOUT - optional output variable to store the value of the gradient
           of F at the minimum.

     If the multivariate function has  more than one minimum, which minimum
     is  returned  is  undefined  (although  it  depends  on  the  starting
     parameters X).

     In  case of  early  termination, the  best  solution found  so far  is
     returned.


   KEYWORDS
     EXTRA - Supplemental  argument  for F;  if  non-nil,  F  is called  as
         F(X,GX,EXTRA) so its prototype must be: func F(x, &gx, extra).
     XMIN, XMAX  - Lower/upper bounds for  X.  Must be  conformable with X.
         For  instance  with  XMIN=0,  the non-negative  solution  will  be
         returned.
     METHOD - Scalar integer which  defines the optimization method to use.
         Conjugate  gradient   algorithm  is  used  if  one   of  the  bits
         OP_FLAG_POLAK_RIBIERE,         OP_FLAG_FLETCHER_REEVES,         or
         OP_FLAG_HESTENES_STIEFEL  is  set;  otherwise,  a  limited  memory
         variable  metric algorithm  (VMLM-B) is  used.  If  METHOD  is not
         specified and  if MEM=0, a conjugate gradient  search is attempted
         with flags: (OP_FLAG_UPDATE_WITH_GP |
                      OP_FLAG_SHANNO_PHUA    |
                      OP_FLAG_MORE_THUENTE   |
                      OP_FLAG_POLAK_RIBIERE  |
                      OP_FLAG_POWELL_RESTART)
         otherwise VMLM-B is used with flags: (OP_FLAG_UPDATE_WITH_GP |
                                               OP_FLAG_SHANNO_PHUA    |
                                               OP_FLAG_MORE_THUENTE).
         See documentation  of op_get_flags to  figure out the  allowed bit
         flags and their meaning.
     MEM -  Number of previous  directions used in variable  metric limited
         memory   method  (default   min(5,  numberof(X))).    If   MEM  is
         explicitely set to zero, the default method is conjugate-gradient.
     MAXITER - Maximum number of iterations (default: no limits).
     MAXEVAL - Maximum number of function evaluations (default: no limits).
     FTOL -  Relative function  change tolerance for  convergence (default:
         1.5e-8).
     GTOL - Gradient tolerance  for  convergence (default: 3.7e-11).
     VERB - Verbose mode?  If non-nil and non-zero, print out  information
        every VERB iterations and for the final one.
     SFTOL, SGTOL, SXTOL, SXBIG - Line   search   tolerance  and  safeguard
        parameters (see op_csrch).

   SEE ALSO: op_get_flags, op_csrch,
             op_cgmnb_setup, op_cgmnb_next,
             op_vmlmb_setup, op_vmlmb_next. */
{
  /* Setup options. */
  use_extra   = ! is_void(extra);
  use_maxiter = ! is_void(maxiter);
  use_maxeval = ! is_void(maxeval);
  cmplx = (structof(x) == complex);
  if (cmplx) {
    if (! is_void(xmin) || ! is_void(xmax)) {
      error, "bounds forbidden for complex type";
    }
  } else {
    if (! is_void(xmin)) x = max(x, xmin);
    if (! is_void(xmax)) x = min(x, xmax);
  }

  /* Initializes workspace according to optimization method and options. */
  if (is_void(mem)) mem = min(5, numberof(x));
  ws = op_mnb_setup(method, mem,
                    gtol=gtol, ftol=ftol,
                    sftol=sftol, sgtol=sgtol, sxtol=sxtol, sxbig=sxbig);
  if (op_get_flags(ws) & OP_ANY_CONJUGATE_GRADIENT_METHOD) {
    next_stage = op_cgmnb_next;
  } else {
    next_stage = op_vmlmb_next;
  }
  stage = 0;

  msg = string();
  if (verb) {
    elapsed = array(double, 3);
    timer, elapsed;
    cpu_start = elapsed(1);
  }
  
  elapsed = array(double, 3);
  timer, elapsed;
  cpu_start = elapsed(1);
  
  eval = 0;
  for (;;) {
    local gx; /* to store the gradient */
    if (stage <= 1) {
      /* Compute function value and gradient at X. */
      fx = (use_extra ? f(x, gx, extra) : f(x, gx));
      ++eval;
    }
    if (stage != 1) {
      /* Decide to continue/stop and/or display information. */
      iter = op_get_iter(ws);
      if ((stop = (stage >= 3))) {
        msg = op_get_msg(ws);
      } else if ((stop = (use_maxiter && iter > maxiter))) {
        msg = swrite(format="warning: too many iterations (%d)", iter);
      } else if ((stop = (use_maxeval && eval > maxeval))) {
        msg = swrite(format="warning: too many function evaluations (%d)",
                     eval);
      }
      if (verb) {
        if (eval == 1) {
          write, format="%s  %s\n%s  %s\n",
            " ITER    EVAL     CPU [s]            FUNC             max(|G|)",
            " STEPLEN",
            "------  ------  ----------  -----------------------  ---------",
            "---------";
        }
        if (stop || iter%verb == 0) {
          timer, elapsed;
          cpu = elapsed(1) - cpu_start;
          write, format="%6d  %6d %11.2e  %23.16e  %9.1e  %9.1e\n",
            iter, eval, cpu, fx, op_get_amaxg(ws), op_get_step(ws);
        }
      }
      if (stop) {
        if (verb && msg) write, format="%s\n", msg;
        return x;
      }
    }

    /* Call optimizer. */
    if (cmplx) {
      x = op_z2d(x);
      gx = op_z2d(gx);
      stage = next_stage(x, fx, gx, ws, xmin=xmin, xmax=xmax);
      x = op_d2z(x);
      gx = op_d2z(gx);
    } else {
      stage = next_stage(x, fx, gx, ws, xmin=xmin, xmax=xmax);
    }
    
    timer, elapsed;
    cpu = elapsed(1) - cpu_start;
    fout=fx;
    gout=gx;
    cpuout=cpu;
    evalout=eval;
    steplenout=op_get_step(ws);
    maxgout=op_get_amaxg(ws);
  }
}

func op_mnb_next(&x, &f, &g, ws, xmin=, xmax=)
{
  return ((op_get_flags(ws) & OP_ANY_CONJUGATE_GRADIENT_METHOD) ?
          op_cgmnb_next : op_vmlmb_next)(x, f, g, ws, xmin=xmin, xmax=xmax);
}

func op_mnb_setup(method, mem, gtol=, ftol=, sftol=, sgtol=, sxtol=, sxbig=)
{
  if (is_void(method)) method = (mem==0 ?
                                 (OP_FLAG_UPDATE_WITH_GP |
                                  /*OP_FLAG_SHANNO_PHUA |*/
                                  OP_FLAG_MORE_THUENTE |
                                  OP_FLAG_POLAK_RIBIERE |
                                  OP_FLAG_POWELL_RESTART):
                                 (OP_FLAG_UPDATE_WITH_GP |
                                  OP_FLAG_SHANNO_PHUA |
                                  OP_FLAG_MORE_THUENTE));
  if ((method & OP_ANY_CONJUGATE_GRADIENT_METHOD) != 0) {
    /* Use conjugate gradients. */
    return op_cgmnb_setup(flags=method, gtol=gtol, ftol=ftol,
                          sftol=sftol, sgtol=sgtol, sxtol=sxtol, sxbig=sxbig);
  }

  /* Use variable metric. */
  return op_vmlmb_setup((is_void(mem) ? 5 : mem),
                        flags=method, gtol=gtol, ftol=ftol,
                        sftol=sftol, sgtol=sgtol, sxtol=sxtol, sxbig=sxbig);
}

/*---------------------------------------------------------------------------*/
/* TESTS AND CHECKING */

func op_check_gradient(f, g, x0, number=, tiny=, dir=)
/* DOCUMENT op_check_gradient(f, g, x)
     Compare gradient function with gradient estimated by finite
     differences.  F is the function, G is the gradient and X the
     parameters.

     The number of gradient values to check may be specified by keyword
     NUMBER (the subset of parameters is randomly chosen). The default is
     to compute the finite difference gradient for all parameters.

     The (relative) finite difference step size cand be specified by
     keyword TINY (default 1e-5).

     Keyword DIR can be used to specify which kind of finite differences to
     use to estimate the gradient: forward (DIR>0), backward (DIR<0) or
     centered (DIR=0 or nil).  The default is to use centered finite
     differences which are more precise (of order TINY^3) but twice more
     expensive to compute.

   SEE ALSO: . */
{
  n = numberof(x0);
  if (is_void(number) || number >= n) {
    list = indgen(n);
    number = n;
  } else {
    list = long(1 + n*random(number));
  }
  if (is_void(tiny)) tiny = 1e-5;
  x1 = x0;
  f0 = f(x0);
  g0 = g(x0)(list);
  g1 = array(double, number);
  if (dir) {
    /* Forward or backward differences. */
    tiny = abs(tiny)*sign(dir);
    for (i=1 ; i<=number ; ++i) {
      l = list(i);
      x = x0(l);
      x1(l) = x + (dx = tiny*max(1.0, abs(x)));
      g1(i) = (f(x1) - f0)/dx;
      x1(l) = x;
    }
  } else {
    /* Centered differences. */
    for (i=1 ; i<=number ; ++i) {
      l = list(i);
      x = x0(l);
      dx = tiny*max(1.0, abs(x));
      x1(l) = x + dx;
      f1 = f(x1);
      x1(l) = x - dx;
      f2 = f(x1);
      x1(l) = x;
      g1(i) = (f1 - f2)/dx; 
    }
    g1 *= 0.5;
  }
  abs_err = abs(g1 - g0);
  max_abs_err = max(abs_err);
  avg_abs_err = avg(abs_err);
  rms_abs_err = sqrt(avg(abs_err*abs_err));
  u = abs(g0) + abs(g1);
  rel_err = 2*(g1 != g2)*abs_err/(u + !u);
  max_rel_err = max(rel_err);
  avg_rel_err = avg(rel_err);
  rms_rel_err = sqrt(avg(rel_err*rel_err));

  write, format="GRADIENT CHECK WITH: tiny=%.1e,  number=%d\n",
    double(tiny), number;
  write, format="ABSOLUTE ERROR: max=%.1e,  avg=%.1e,  rms=%.1e\n",
    max_abs_err, avg_abs_err, rms_abs_err;
  write, format="RELATIVE ERROR: max=%.1e,  avg=%.1e,  rms=%.1e\n",
    max_rel_err, avg_rel_err, rms_rel_err;
}

/*---------------------------------------------------------------------------*/
/* MAIN DRIVERS FOR LIMITED MEMORY MINIMIZATION */

/*
 * The workspace is an array of pointers (memory needed):
 *   WS(1) = &ISAVE
 *   WS(2) = &DSAVE
 *   WS(3) = &MSG
 *   WS(4) = &D       anti-search direction               (N)
 *   WS(5) = &XBEST   best solution found so far          (0/N)
 *   WS(6) = &GBEST   gradient at XBEST                   (0/N)
 *   WS(7:) = specific to minimization method
 *
 *   ISAVE(1) = BRACKT     bracket in line search?
 *   ISAVE(2) = LS_STAGE   stage in line search
 *   ISAVE(3) = STAGE      stage in multi-variate optimization
 *   ISAVE(4) = ITER       number of iterations
 *   ISAVE(5) = RESTARTS   number of restarts with steepest descent
 *   ISAVE(6) = REJECTS    number of update rejections
 *   ISAVE(7) = FLAGS      bitwise flags
 *
 *   DSAVE(1:13) is used by op_csrch
 *   DSAVE(14) = SFTOL
 *   DSAVE(15) = SGTOL
 *   DSAVE(16) = SXTOL
 *   DSAVE(17) = SXBIG
 *   DSAVE(18) = FTOL
 *   DSAVE(19) = GTOL
 *   DSAVE(20) = F0
 *   DSAVE(21) = DG
 *   DSAVE(22) = DG0
 *   DSAVE(23) = STP
 *   DSAVE(24) = STPMIN
 *   DSAVE(25) = STPMAX
 *   DSAVE(26) = AMAXG      maximum absolute value of (projected) gradient
 *   DSAVE(27) = FBEST      function value at best solution found so far
 *   DSAVE(28) = STPBEST    step length at best solution found so far
 *   DSAVE(29) = AMAXGBEST  step length at best solution found so far
 *
 * FIXME: it is possible to save some memory (GBEST/GINIT can be the same?)
 *
 * To check worspace numbers:
 *   grep ws\( yorick/optimpack.i|sed 's/^.*\(ws([0-9])\).*$/\1/'|sort -u
 */

extern op_get_msg; // link to builtin function
func op_get_msg(ws) { return *ws(3); }
/* DOCUMENT op_get_msg(ws)
     Get message stored into workspace WS.
   SEE ALSO: op_cgmnb_setup, op_vmlmb_setup. */

local OP_STAGE_START, OP_STAGE_FG,   OP_STAGE_NEWX;
local OP_STAGE_CONV,  OP_STAGE_WARN, OP_STAGE_ERROR;
extern op_get_stage;
func op_get_stage(ws) { return (*ws(1))(3); }
/* DOCUMENT op_get_stage(ws)
     Get stage (from workspace WS) in multidimensional minimization:
        0 (OP_STAGE_START) - Initialization: caller must choose a starting
           solution and evaluate the function and gradient.
        1 (OP_STAGE_FG) - Caller must evaluate the function and gradient.
        2 (OP_STAGE_NEWX) - New iterate has been computed, improved solution,
           function and gradient are available for examination.
        3 (OP_STAGE_CONV) - Algorithm has converged.  The symbols storing
           the variables, the function value and the gradient have been set
           to contain the solution.
        4 (OP_STAGE_WARN) - Warning: the algorithm is not able to satisfy
           the convergence conditions.  The symbols storing the variables,
           the function value and the gradient have been set to contain the
           best approximation found so far.
        5 (OP_STAGE_ERROR) - There is an error in  the input arguments.
     If STAGE >= 3, you can use op_get_msg (which see) to retrieve the
     corresponding message.

   SEE ALSO: op_cgmnb_setup, op_vmlmb_setup. */
OP_STAGE_START = 0; /* first entry, start search */
OP_STAGE_FG    = 1; /* computation of F and G requested */
OP_STAGE_NEWX  = 2; /* new improved solution available for inspection */
OP_STAGE_CONV  = 3; /* search has converged */
OP_STAGE_WARN  = 4; /* search aborted with warning */
OP_STAGE_ERROR = 5; /* search aborted with error */

/* "local" declarations must fit in the 10 lines above DOCUMENT comment */
extern op_get_flags;
func op_get_flags(ws) { return (*ws(1))(7); }
local OP_FLAG_ARMIJO, OP_FLAG_MORE_THUENTE, OP_FLAG_SHANNO_PHUA;
local OP_FLAG_FLETCHER_REEVES, OP_FLAG_POLAK_RIBIERE;
local OP_FLAG_HESTENES_STIEFEL, OP_FLAG_POWELL_RESTART;
local OP_ANY_CONJUGATE_GRADIENT_METHOD, OP_FLAG_LBFGS;
local OP_FLAG_SOLVE_WITH_GP;
local OP_FLAG_RESTRICT_SOLVE;
local OP_FLAG_UPDATE_WITH_GP;
/* DOCUMENT op_get_flags(ws)

     Get value of bitwise  flags in multidimensional minimization workspace
     WS.  The allowed flags are:

     OP_FLAG_ARMIJO - Use Armijo's line search (only check for sufficient
         descent condition).

     OP_FLAG_MORE_THUENTE - Use Mori & Thuente's line search (sufficient
	 decrease and curvature conditions).

     OP_FLAG_SHANNO_PHUA - Use Shanno & Phua's formula to pre-scale the
	 search direction.

     OP_FLAG_FLETCHER_REEVES - Use Fletcher & Reeves's formula in
         conjugate-gradient recursion (ok if the objective function is
         quadratic, for non-linear function, Polak & Ribihre is superior).

     OP_FLAG_POLAK_RIBIERE - Use Polak & Ribihre's formula in
         conjugate-gradient recursion (generaly better than Fletcher &
         Reeves for non-quadratic objective function).

     OP_FLAG_HESTENES_STIEFEL - Use Hestenes & Stiefel 's formula in
         conjugate-gradient recursion.

     OP_FLAG_POWELL_RESTART - Use Powell's rule to restart the
         conjugate-gradient recursion whenever the parameter BETA becomes
         non-positive.

     OP_FLAG_LBFGS - Use limited memory Broyden-Fletcher-Goldfarb-Shanno
         (BFGS) variable metric method.

     OP_FLAG_SOLVE_WITH_GP - Apply approximation of inverse Hessian to the
	 projected gradient rather than to the gradient.

     OP_FLAG_RESTRICT_SOLVE - Recursive formula to approximate the inverse
	 Hessian is applied to the set of free parameters only.  In this
	 case, there is no need to specify OP_FLAG_SOLVE_WITH_GP.

     OP_FLAG_UPDATE_WITH_GP - Update Hessian approximation with the
         projected gradients rather than with the gradients.

     OP_FLAG_SAVE_BEST - Save best solution found so far in workspace.
         This may (or may not!) improve convergence to the cost of some
         extra memory.

     Global variable OP_ANY_CONJUGATE_GRADIENT_METHOD is set with bits
     OP_FLAG_POLAK_RIBIERE, OP_FLAG_FLETCHER_REEVES and
     OP_FLAG_HESTENES_STIEFEL.

   SEE ALSO: op_cgmnb_setup, op_multi, op_vmlmb_setup. */

/* Following values must match those in optimpack.h */
OP_FLAG_ARMIJO           = (1<< 0);
OP_FLAG_MORE_THUENTE     = (1<< 1);
OP_FLAG_SHANNO_PHUA      = (1<< 2);
OP_FLAG_SOLVE_WITH_GP    = (1<< 3);
OP_FLAG_RESTRICT_SOLVE   = (1<< 4);
OP_FLAG_UPDATE_WITH_GP   = (1<< 5);
OP_FLAG_FLETCHER_REEVES  = (1<< 6);
OP_FLAG_POLAK_RIBIERE    = (1<< 7);
OP_FLAG_HESTENES_STIEFEL = (1<< 8);
OP_FLAG_POWELL_RESTART   = (1<< 9);
OP_FLAG_SAVE_BEST        = (1<<10);
OP_ANY_CONJUGATE_GRADIENT_METHOD = (OP_FLAG_POLAK_RIBIERE |
                                    OP_FLAG_FLETCHER_REEVES |
                                    OP_FLAG_HESTENES_STIEFEL);

extern op_get_iter;
func op_get_iter(ws) { return (*ws(1))(4); }
/* DOCUMENT op_get_iter(ws)
     Get number of iterations (i.e. number of successul improvements of the
     solution).  WS is the workspace object used in multidimensional
     minimization.
   SEE ALSO: op_cgmnb_setup, op_vmlmb_setup. */

extern op_get_restarts;
func op_get_restarts(ws) { return (*ws(1))(5); }
/* DOCUMENT op_get_restarts(ws)
     Get number of restarts with the steepest descent.  WS is the workspace
     object used in multidimensional minimization.
   SEE ALSO: op_cgmnb_setup, op_vmlmb_setup. */

extern op_get_rejects;
func op_get_rejects(ws) { return (*ws(1))(6); }
/* DOCUMENT op_get_rejects(ws)
     Get number of updates of Hessian approximation which have been
     rejected.  WS is the workspace object used in multidimensional
     minimization.
   SEE ALSO: op_cgmnb_setup, op_vmlmb_setup. */

extern op_get_step;
func op_get_step(ws) { return (*ws(2))(23); }
/* DOCUMENT op_get_step(ws)
     Get step length in line search.  WS is the workspace object used in
     multidimensional minimization.
   SEE ALSO: op_cgmnb_setup, op_vmlmb_setup. */

extern op_get_amaxg;
func op_get_amaxg(ws) { return (*ws(2))(26); }
/* DOCUMENT op_get_amaxg(ws)
     Get maximum absolute value of (projected) gradient.  WS is the
     workspace object used in multidimensional minimization.

   SEE ALSO: op_cgmnb_setup, op_vmlmb_setup. */

func _op_setup(nws, nisave, ndsave)
/* DOCUMENT _op_setup(nws, nisave, ndsave)
     Private helper routine to prepare workspaces for limited-memory
     optimization algorithms.  Arguments are:
        NWS    - number of slots in worspace (NWS >= 6)
        NISAVE - number of integer values (NISAVE >= 7)
        NDSAVE - number of real values (NDSAVE >= 29)
     External symbols (which must be "local" for the caller):
       _op_message                - error message (when nil is returned)
       flags                      - bitwise flags
       ftol, gtol                 - global convergence criteria
       sftol, sgtol, sxtol, sxbig - line search parameters

     Default values: FTOL = EPSILON^(1/2) = 1.5e-8
                     GTOL = EPSILON^(2/3) = 3.7e-11

   SEE ALSO: _op_next, op_cgmnb_setup, op_vmlmb_setup. */
{
  {extern _op_message;}
  {extern flags, ftol, gtol, sftol, sgtol, sxtol, sxbig;}
  not_scalar_real    = op_not_scalar_real;
  not_scalar_integer = op_not_scalar_integer;

  if (not_scalar_real(ftol, 1.5e-8) || ftol < 0.0) {
    _op_message = "bad value for keyword FTOL";
  } else if (not_scalar_real(gtol, 3.7e-11) || gtol < 0.0) {
    _op_message = "bad value for keyword GTOL";
  } else if (not_scalar_real(sftol, 0.001) || sftol <= 0.0 || sftol >= 1.0) {
    _op_message = "bad value for keyword SFTOL";
  } else if (not_scalar_real(sgtol, 0.9) || sgtol <= 0.0 || sgtol >= 1.0) {
    _op_message = "bad value for keyword SGTOL";
  } else if (sftol >= sgtol) {
    _op_message = swrite(format="SFTOL (%g) must be less than SGTOL (%g)",
                            sftol, sgtol);
  } else if (not_scalar_real(sxtol, 0.1) || sxtol <= 0.0 || sxtol >= 1.0) {
    _op_message = "bad value for keyword SXTOL";
  } else if (not_scalar_real(sxbig, 1e40) || sxbig < 2.0) {
    _op_message = "bad value for keyword SXBIG";
  } else {
    /* Setup workspace. */
    isave = array(long, nisave);
    isave(5) = -1; /* number of restarts */
    isave(7) = flags;

    dsave = array(double, ndsave);
    dsave(14) = sftol;
    dsave(15) = sgtol;
    dsave(16) = sxtol;
    dsave(17) = sxbig;
    dsave(18) = ftol;
    dsave(19) = gtol;
    dsave(26) = -1; /* AMAXG */

    ws = array(pointer, nws);
    ws(1) = &isave;
    ws(2) = &dsave;
    return ws;
  }
}

func _op_next(xinit, ginit, next_dir, save_state)
/* DOCUMENT _op_next(xinit, ginit, next_dir, save_state)
     Private driver.
   SEE ALSO: op_cgmnb_next, op_vmlmb_next. */
{
  {extern f, x, g, ws, xmin, xmax;}
  local isave; eq_nocopy, isave, *ws(1);
  local dsave; eq_nocopy, dsave, *ws(2);
  local d;     eq_nocopy, d,     *ws(4);

  /* Provide an interpreted version of the unref function if not
     builtin in current interpreter. */
  if (is_func(unref) != 2) unref = op_unref;

  /* Restore local variables. */
  stage  = isave(3);
  flags  = isave(7);
  sftol  = dsave(14);
  sgtol  = dsave(15);
  sxtol  = dsave(16);
  sxbig  = dsave(17);
  ftol   = dsave(18);
  gtol   = dsave(19);
  f0     = dsave(20);
  dg     = dsave(21);
  dg0    = dsave(22);
  stp    = dsave(23);
  stpmin = dsave(24);
  stpmax = dsave(25);
  amaxg  = dsave(26);

  /* Figure out the set of free variables (FIXME: this computation can be
     avoided for some stages if array FREE is saved in workspace WS). */
  if (is_void(xmin)) {
    if (is_void(xmax)) {
      /* no bounds */
      constrained = 0n;
    } else {
      /* only upper bounds */
      constrained = 1n;
      free = ((x < xmax) | (g > 0.0));
    }
  } else {
    if (is_void(xmax)) {
      /* only lower bounds */
      constrained = 1n;
      free = ((x > xmin) | (g < 0.0));
    } else {
      /* lower and upper bounds */
      constrained = 1n;
      free = (((x > xmin) | (g < 0.0)) & ((x < xmax) | (g > 0.0)));
    }
  }
  if (constrained && noneof(free)) {
    f0 = f;
    dg0 = dg = 0.0;
    msg = "convergence: no free variables";
    stage = 3;
    goto done; 
  }

  if (stage <= 1) {
    /* Compute maximum magnitude of (projected) gradient. */
    if (stage || amaxg < 0) amaxg = max(abs((constrained ? free*g : g)));

    /* Maybe save best solution found so far. */
    if ((flags & OP_FLAG_SAVE_BEST) &&
        (stage == 1 ? f < dsave(27) : isave(5) /* restarts */ == -1)) {
      copy = x; ws(5) = &copy; /* xbest */
      copy = g; ws(6) = &copy; /* gbest */
      dsave(27) = f; /* fbest */
      dsave(28) = (stage == 1 ? stp : 0.0); /* stpbest */
      dsave(29) = amaxg; /* amaxgbest */
    }
  }

  if (stage == 1) {
    /* Line search in progress. */
    dg = -sum(d*(constrained ? free*g : g));
    task = "fg"; /* continue with line search */
  } else {
    /* Compute new search direction. */
    if (stage == 2 || stage == 3) {
      /* Update Hessian information and use recursion to guess next
         search direction. */
      ws(4) = pointer(); /* free memory, previous direction still in D */
      stage = next_dir();
      if (stage == 1) {
        /* Memorize new search direction. */
        ws(4) = &d;
      } else if (stage == 0) {
        /* Recursion failed to provide a descent direction,
           restart with steepest descent. */
        d = [];
      } else {
        goto done;
      }
    } else if (stage) {
      /* Make sure the algorithm is re-started if caller insists on
         continuing after an error or a warning. */
      stage = 0;
    }
    if (! stage) {
      /* Start/restart the algorithm.  Use steepest descent. */
      ++isave(5); /* restarts */
      d = double((constrained ? g*free : g)); /* FIXME: pre-conditioning */
      if ((q = op_nrm2(d)) <= 0.0) {
        f0 = f;
        dg0 = dg = 0.0;
        msg = swrite(format= "convergence: local minimum found %s",
                     (constrained ?
                      "which exactly satisfies K|hn-Tucker conditions" :
                      "whith norm of gradient equal to zero"));
        stage = 3;
        goto done; 
      }
      d = (1.0/q)*unref(d);
      ws(4) = &d;
      dg = -q;
    }

    /* Setup for start of line search. */
    save_state;
    f0 = f;
    dg0 = dg;
    stpmax = op_get_max_step(-sxbig, *ws(4), xmin, x, xmax);
    if (stpmax < 0.0) {
      msg = "error: some parameters violate their bounds";
      stage = 5;
      goto done;
    }
    stp = min(1.0, stpmax);
    task = "start"; /* start line search */
  }

  /* Call line search iterator. */
  op_csrch,f,dg,stp,sftol,sgtol,sxtol,stpmin,stpmax,task,isave,dsave;  
  if (task == "fg") {
    /* Compute the new iterate.  Next call will continue line search. */
    x = op_apply_bounds(xmin, xinit(ws) - stp*d, xmax);
    msg = string();
    stage = 1;
  } else if ((ident = strpart(task, 1:4) == "conv") ||
             task == "warning: sxtol test satisfied") {
    /* Line search has converged (or relative size of interval of
       uncertainty is small enough).  Test for overall convergence
       otherwise set STAGE to signal a new iterate.  In any case,
       next call will compute a new L-BFGS search direction. */
    if (amaxg <= (q = gtol*(1.0 + abs(f)))) {
      /* K|hn-Tucker conditions satisfied within tolerance. */
      msg = swrite(format="convergence: maximum absolute value of %sgradient less or equal %g", (constrained?"projected ":""), q);
      stage = 3;
    } else if (max(abs(f - f0), abs(stp*dg0)) <= ftol*max(abs(f0), abs(f))) {
      msg = swrite(format="convergence: relative function reduction less or equal FTOL=%g", ftol);
      stage = 3;
    } else {
      msg = "new estimate available for inspection";
      stage = 2;
    }
    ++isave(4); /* iter */
    if ((flags & OP_FLAG_SAVE_BEST) && (dsave(27) /* fbest */ < f)) {
      /* Restore best solution found so far. */
      x = *ws(5); /* xbest */
      g = *ws(6); /* gbest */
      f = dsave(27); /* fbest */
      stp = dsave(28); /* stpbest */
      amaxg = dsave(29); /* amaxgbest */
    }
  } else {
    /* Error or warning.  Algorithm will be completely restarted if caller
       insists on continuing. */
    msg = task;
    stage = (ident == "warn" ? 4 : 5);
    if (f < f0) {
      /* Current solution is better than previous one. */
      ++isave(4); /* iter */
    } else {
      /* Restore better solution (note that with update_with_gp, the
         restored gradient is the projected one but this probably do not
         matter). */
      if (flags & OP_FLAG_SAVE_BEST) {
        x = *ws(5); /* xbest */
        g = *ws(6); /* gbest */
        f = dsave(27); /* fbest */
        stp = dsave(28); /* stpbest */
        amaxg = dsave(29); /* amaxgbest */
      } else {
        x = xinit(ws);
        g = ginit(ws);
        f = f0;
        amaxg = max(abs(g)); // FIXME: MUST BE PROJECTED GRADIENT HERE
      }
    }
  }

  /* Save local variables (but the constant ones). */
 done:
  ws(3) = &msg;
  isave(3) = stage;

  dsave(20) = f0;
  dsave(21) = dg;
  dsave(22) = dg0;
  dsave(23) = stp;
  dsave(24) = stpmin;
  dsave(25) = stpmax;
  dsave(26) = amaxg;
  return stage;
}

/*---------------------------------------------------------------------------*/
/* CONJUGATE GRADIENT MINIMIZATION */

/*     Gilbert, J.C. and Nocedal, J. (1992). "Global Convergence Properties 
       of Conjugate Gradient Methods", SIAM Journal on Optimization, Vol. 2,
       pp. 21-42.


       D(k) is the anti-search direction:

          x(k+1) = x(k) - alpha(k)*d(k)

       with (conjugate-gradient)

          d(k) = g(k) + beta(k)*d(k-1)

       Shanno-Phua's formula for length of trial step along the new
       search direction:

          alpha(k)_trial = alpha(k) d(k-1)'.g(k-1) / d(k)'.g(k)

       means that: alpha(k) x d(k)'.g(k)  remains almost constant


*/

/* Conjugate-gradient minimization:
 *   WS(1:6) see _op_setup
 *   WS(7) = &X0 parameters at start of line search
 *   WS(8) = &G0 (projected) gradient at X0
 */
extern op_cgmnb_setup;
func op_cgmnb_setup(nil, flags=, ftol=, gtol=,
                    sftol=, sgtol=, sxtol=, sxbig=)
{
  local _op_message;
  if (is_void(flags)) flags = (OP_FLAG_POLAK_RIBIERE|
                               OP_FLAG_POWELL_RESTART);
  if (is_void(sftol)) sftol = 1e-4;
  if (is_void(sgtol)) sgtol = 1e-1; /* between 0.1 and 0.01 */
  if (is_void(sxtol)) sxtol = 1e-17;
  ws = _op_setup(8, 7, 29);
  if (is_void(ws)) error, _op_message;
  return ws;
}

extern op_cgmnb_next;
func op_cgmnb_next(&x, &f, &g, ws, xmin=, xmax=)
{
  return _op_next(_op_cgmnb_xinit,
                  _op_cgmnb_ginit,
                  _op_cgmnb_next_dir,
                  _op_cgmnb_save_state);
}

func _op_cgmnb_xinit(ws) { return *ws(7); }
func _op_cgmnb_ginit(ws) { return *ws(8); }
func _op_cgmnb_save_state {
  {extern x, f, g, ws, isave, dsave;}
  ws(7) = ws(8) = pointer(); /* free some memory */
  x0 = double(x); ws(7) = &x0; /* store private copy of X */
  g0 = double(g); ws(8) = &g0; /* store private copy of G */
}
/* next_dir(d) returns next (anti) search direction given the previous one D
   returns nil if restart is needed */
func _op_cgmnb_next_dir(nil)
{
  {extern x, f, g, ws, isave, dsave, d, dg, stp, dg0;}
  local g0; eq_nocopy, g0, _op_cgmnb_ginit(ws); /* previous gradient */
  if (flags & OP_FLAG_FLETCHER_REEVES) {
    beta = sum(g*g)/sum(g0*g0);
  } else if (flags & OP_FLAG_POLAK_RIBIERE) {
    beta = sum((g - g0)*g)/sum(g0*g0);
  } else if (flags & OP_FLAG_HESTENES_STIEFEL) {
    tmp = g - g0;
    beta = -sum(tmp*g)/sum(tmp*d);
    tmp = [];
  } else {
    /* use Polak & Ribihre's formula by default. */
    beta = sum((g - g0)*g)/sum(g0*g0);
  }
  if (((flags & OP_FLAG_POWELL_RESTART) ? beta <= 0.0 : beta == 0.0)) {
    /* Invalid BETA: restart algorithm. */
    return 0;
  }
  d = g + beta*unref(d);
  if ((dg = -sum(d*g)) >= 0.0) {
     /* Recursion failed to provide a descent direction: restart algorithm. */
    return 0;
  }
  if (flags & OP_FLAG_SHANNO_PHUA) {
    /* Use Shanno-Phua's formula for length of trial step along the new
       search direction. */
    // FIXME: does not work at all...
    gamma = stp*dg0/dg;
    //write, format="gamma=%g; beta=%g;\n", gamma, beta;
    if (gamma != 1) {
      d = gamma*unref(d);
      dg *= gamma;
    }
  }
  return 1;
}

/*---------------------------------------------------------------------------*/
/* LIMITED MEMORY VARIABLE METRIC MINIMIZATION */

/* L-LBFGS minimization:
 *   WS(1:6) see _op_setup
 *   WS(7) = &S   parameters differences
 *   WS(8) = &Y   (projected) gradient differences
 *   WS(9) = &RHO
 *
 * ISAVE is an array of 6 long's:
 *   ISAVE(1:2) is used by op_csrch
 *   ISAVE(3:7) is used by _op_next
 *   ISAVE( 8) = M    maximum number of memorized update pairs
 *   ISAVE( 9) = MP   effective number of memorized update pairs
 *   ISAVE(10) = MARK location of last memorized update pair
 */

extern op_vmlmb_setup;
func op_vmlmb_setup(m, flags=, ftol=, gtol=,
                    sftol=, sgtol=, sxtol=, sxbig=)
/* DOCUMENT ws = op_vmlmb_setup(m)

     Returns a new workspace for multidimensional iterative minimization
     by limited memory variable metric algorithm VMLMB.  M is the number of
     updates to memorize. FIXME:

   SEE ALSO:
 */
{
  local _op_message;
  if (is_void(flags)) flags = (OP_FLAG_UPDATE_WITH_GP|
                               OP_FLAG_SHANNO_PHUA|
                               OP_FLAG_MORE_THUENTE);
  if (is_void(sftol)) sftol = 0.001;
  if (is_void(sgtol)) sgtol = 0.9;
  if (is_void(sxtol)) sxtol = 0.1;
  if (op_not_scalar_integer(m) || m < 1) error, "bad value for M";
  ws = _op_setup(9, 10, 29);
  if (is_void(ws)) error, _op_message;
  ws(7) = &array(pointer, m);  // S
  ws(8) = &array(pointer, m);  // Y
  ws(9) = &array(double, m);   // RHO
  (*ws(1))(8) = m; /* ISAVE(8) = M */
  (*ws(1))(10) = -1; /* ISAVE(10) = MARK */
  return ws;
}

extern op_vmlmb_next;
func op_vmlmb_next(&x, &f, &g, ws, xmin=, xmax=, free=)
/* DOCUMENT stage = op_vmlmb_next(x, fx, gx, ws);

     Performs next  step in limited-memory variable  metric method possibly
     with  simple  bound constraints  (VMLM-B).   VMLM-B  computes a  local
     minimizer of  a function of N  variables by a  limited memory variable
     metric (BFGS) method; optionally,  the parameters may be bounded.  The
     user must evaluate the function and the gradient.

     VMLM-B   is  implemented   via  two   functions:   op_vmlmb_setup  for
     initialization  and  op_vmlmb_next   for  further  iterations.   These
     functions use reverse communication.   The user must choose an initial
     approximation  X  to the  minimizer,  evaluate  the  function and  the
     gradient at X, and make the  initial call with a fresh workspace given
     by  op_vmlmb_setup.  On  exit STAGE  indicates the  required  action -
     alternatively,  STAGE  can   be  obtained  by  op_get_stage(WS).   See
     op_get_stage for the meaning of the returned value.

     X  is a  double  precision  array of  length  N.  On  entry,  X is  an
     approximation to the solution.  On  exit with STAGE = OP_STAGE_CONV, X
     is the final approximation.

     FX is the value  of the function at X.  On final  exit, FX is set with
     the function value at the final approximation X.

     GX is a double precision array of length N.  On entry, GX is the value
     of the gradient at X.  On final  exit, GX is the value of the gradient
     at the final approximation X.

     WS is a workspace array  as returned by op_vmlmb_setup (which see) and
     must not be changed between calls to op_vmlmb_next.

     Keywords XMIN and/or  XMAX can be set to specify  a lower and/or upper
     bounds for the parameters.  If not  nil, the value of XMIN and/or XMAX
     must be a scalar (same bound  for all parameters) or an array with the
     same  number of elements  as X  (elementwise bounds).   If there  is a
     mixture of bounded and unbounded  parameters, the value of elements of
     XMIN and/or  XMAX corresponding to unbounded parameters  should be set
     to -HUGE or +HUGE respectively where  HUGE is a very large value which
     does not overflow (e.g.  HUGE=1.79769e+308 for IEEE double precision).
     Note that  if the caller uses  and properly instanciates  the value of
     keyword  FREE, XMIN and  XMAX are  not required  to account  for bound
     constraints.

     Keyword FREE  can be  used to specify  an optional integer  array with
     length N provided by the caller if the values in X has bounds.  If the
     parameters  have   no  bounds,  FREE  should   be  nil  (unconstrained
     minimization).  Otherwise, elements set  to zero in FREE indicate that
     the corresponding  values in X has  reached a bound and  should not be
     changed during the  next step because the gradient  has the wrong sign
     (i.e.   the  steepest  descent   direction  would  violate  the  bound
     constraints):

             FREE[i] = 0 if i-th value has a lower bound XMIN[i]
                             and X[i]=XMIN[i] and G[i]>=0
                         0 if i-th value has an upper bound XMAX[i]
                             and X[i]=XMAX[i] and G[i]<=0
                         1 (or any non-zero value) otherwise

         FREE  needs only  to be  computed (and  specified) the  first time
         op_vmlmb_next is called and  when STAGE=OP_STAGE_NEWX (i.e.  after a
         successful   step).     FREE   may   also    be   specified   when
         STAGE=OP_STAGE_CONV  (i.e.   after  convergence  if caller  wish  to
         continue with  minimization).  If X has (some)  bounds, the caller
         is responsible for applying the  bounds to X before evaluating the
         function value  F and the gradient G  (i.e. when STAGE=OP_STAGE_FG),
         e.g.:
             if (X[i] < XMIN[i]) X[i] = XMIN[i];
             if (X[i] > XMAX[i]) X[i] = XMAX[i];

         If Q is  not specified (i.e. Q is  NULL) or if Q[i] > 0  for all i
         such that FREE[i] is non-zero, then FREE is left unchanged.

     Q is an optional double precision  array with length N provided by the
         caller and such that diag(Q) is an approximation of the inverse of
         the Hessian matrix.  If Q is NULL, then the inverse of the Hessian
         is approximated by a simple rescaling using Shanno & Phua formula.
         Otherwise, if  FREE is  NULL, all elements  of Q must  be strictly
         greater than zero; else FREE[i] is  set to zero if Q[i] <= 0 (this
         is the  only case where FREE  is modified).  As for  FREE, Q needs
         only to  be specifed  the first time  op_vmlmb is called  and when
         STAGE=2.

   The caller must  not modify the workspace arrays  CSAVE, ISAVE and DSAVE
   between calls to op_vmlmb_setup and further calls to op_vmlmb_next.

   A  typical invocation of  VMLMB for  unconstrained minimization  has the
   following outline:

      // Choose a starting vector:
      for (i=0 ; i<n ; ++i) x[i] = ... ;

      // Allocate and setup workspaces:
      csave = malloc(OP_VMLMB_CSAVE_NUMBER*sizeof(char));
      isave = malloc(OP_VMLMB_ISAVE_NUMBER*sizeof(op_integer_t));
      dsave = malloc(OP_VMLMB_DSAVE_NUMBER(n, m)*sizeof(double));

      ws = op_vmlmb_setup(n, m, gtol, ftol, sbig, sftol, sgtol, sxtol,
                          csave, isave, dsave); // FIXME:
      for (;;) {
        if (stage == OP_STAGE_FG) {
          f = ...;  // evaluate the function at X; store in F
          g = ...;  // evaluate the gradient of F at X; store in G
        } else if (stage == OP_STAGE_NEWX) {
           // New successful step: the approximation X, function F, and
           // gradient G, are available for inspection.
        } else {
          // Convergence, or error, or warning
          fprintf(stderr, "%s\n", csave);
          break;
        }
        // Computes next step:
        stage = op_vmlmb_next(x, &f, g, NULL, NULL, csave, isave, dsave);
      }

   A typical invocation of VMLMB for bound-constrained minimization has the
   following outline:

      // Choose a starting vector:
      for (i=0 ; i<n ; ++i) x[i] = ... ;

      // Allocate and setup workspaces:
      csave = malloc(OP_VMLMB_CSAVE_NUMBER*sizeof(char));
      isave = malloc(OP_VMLMB_ISAVE_NUMBER*sizeof(op_integer_t));
      dsave = malloc(OP_VMLMB_DSAVE_NUMBER(n, m)*sizeof(double));
      stage = op_vmlmb_setup(n, m, gtol, ftol, sbig, sftol, sgtol, sxtol,
                            csave, isave, dsave);
      eval = 0; // number of evaluations
      for (;;) {
        if (stage == OP_STAGE_FG) {

FIXME:  op_bounds_apply

          op_bounds_apply(n, x, xmin, xmax); // aply bound constraints
          f = ...;  // evaluate the function at X; store in F
          g = ...;  // evaluate the gradient of F at X; store in G
          ++eval;
        } else if (stage == OP_STAGE_NEWX) {
           // New successful step: the approximation X, function F, and
           // gradient G, are available for inspection.
        } else {
          // Convergence, or error, or warning
          fprintf(stderr, "%s\n", csave);
          break;
        }
        // Computes next step:
        if (eval == 1 || stage == OP_STAGE_NEWX) {
          // Computes set of free parameters:
          op_bounds_free(n, free, x, g, xmin, xmax);
        }
        stage = op_vmlmb_next(x, &f, g, free, NULL, csave, isave, dsave);
      }


   HISTORY:
     MINPACK-2 Project. April 1995.
     Argonne National Laboratory and University of Minnesota.
     Brett M. Averick, Richard G. Carter, and Jorge J. Mori.

     C-version and improvements (bound constraints, preconditioning, ...).
     February 2003 - March 2003.
     Observatoire de Lyon.
     Eric Thiebaut.
 */
{
  return _op_next(_op_vmlmb_xinit,
                  _op_vmlmb_ginit,
                  _op_vmlmb_next_dir,
                  _op_vmlmb_save_state);
}

func _op_vmlmb_xinit(ws) { return *(*ws(7))((*ws(1))(10)); /* S(MARK) */ }
func _op_vmlmb_ginit(ws) { return *(*ws(8))((*ws(1))(10)); /* Y(MARK) */ }

func _op_vmlmb_save_state
{
  {extern x, f, g, ws, isave, dsave, d, dg, stp, dg0,
     msg, flags, constrained, free;}
  local s; eq_nocopy, s, *ws(7);
  local y; eq_nocopy, y, *ws(8);
  /* Save state at beginning of line search.  (must make _private_ copy
     of X and (projected) G, hence the use of local variables X0 and G0) */
  m    = isave( 8);
  mp   = isave( 9);
  mark = isave(10);
  if (++mark >= m) mark = 0;
  isave(10) = mark;
  s(mark) = y(mark) = pointer(); /* free some memory */
  x0 = double(x);
  g0 = double((constrained && (flags&OP_FLAG_UPDATE_WITH_GP) ? free*g : g));
  s(mark) = &x0;
  y(mark) = &g0;
}

/* next_dir(d) returns next (anti) search direction given the previous one D
   returns nil if restart is needed */
func _op_vmlmb_next_dir(nil)
{
  {extern x, f, g, ws, isave, dsave, d, dg, stp, dg0,
     msg, flags, constrained, free;}
  local s;   eq_nocopy, s,   *ws(7);
  local y;   eq_nocopy, y,   *ws(8);
  local rho; eq_nocopy, rho, *ws(9);
  m    = isave( 8);
  mp   = isave( 9);
  mark = isave(10);
  epsilon = 2.3e-16; /* machine relative precision */

  /* At this stage:
   *   S(mark) = X0 parameters at start of line search
   *   Y(mark) = G0 (projected) gradient at X0
   * we store the effective changes in parameters and (projected)
   * gradient into S(mark) and Y(mark); if the update is rejected, then
   * there will be at most M-1 pairs in the approximation (but this
   * save memory).
   *
   *   MP     is the actual number of memorized pairs
   *   MARK   is the index of the latest pair which is wrapped in
   *          the range [0, ..., M-1] 
   *   MARK-1 is the index of the next latest pair...
   * Note that the same indexing for MARK as in C language is applied; this
   * make the code more similar to the C version and does not matter because
   * index 0 is valid in Yorick (means last one) as far as the wrapping and
   * increment of MARK is correctly done.
   *
   *
   * Updating means:
   *   save (s,y) pair _and_ value of gamma
   *
   * Then we compute Shanno-Phua scale.
   *
   * We apply the 2-loop recursion to compute the search direction.  If
   * the search direction appears to be impossible to compute (not
   * enough valid memorized pairs) or to not be a descent direction, we
   * will use the steepest direction instead.
   *
   * We replace the oldest (s,y) pair with the new starting point and
   * the (projected) gradient.
   */

  /*******************************************************
   ***                                                 ***
   ***  Update L-BFGS approximation of Hessian matrix  ***
   ***                                                 ***
   *******************************************************/

  /* Compute the effective parameter change and (projected) gradient
     change. */
  //write, format="update: mark=%d", mark;
  s(mark) = &((*s(mark)) - x);
  if (constrained && (flags & OP_FLAG_UPDATE_WITH_GP)) {
    /* Compute projected gradient change (projection of the current
       gradient was computed by previous call to this routine and saved
       into array D). */
    y(mark) = &(*y(mark) - free*g);
  } else {
    /* Compute gradient change. */
    y(mark) = &(*y(mark) - g);
  }

  /* Decide to keep or skip the update. */
  yy = sum((*y(mark))*(*y(mark)));
  sy = sum((*s(mark))*(*y(mark)));
  if (yy > 0 && sy > epsilon*yy) {
    /* Accept the L-BFGS update and compute Shanno & Phua's scale factor for
       the last saved pair. */
    if (mp < m) ++mp; /* one more pair */
    rho(mark) = sy;
    gamma = sy/yy; /* Shanno & Phua's scale factor */
  } else {
    /* Skip the L-BFGS update. */
    if (yy == 0 || (sy == 0 && noneof(*s(mark)))) {
      /* No gradient or parameter change: setup workspaces so that we
         will re-start with the speepest descent if caller insists on
         continuing. */
      isave(9) = 0; /* mp */
      isave(10) = -1; /* mark */
      msg = swrite(format="convergence: no %s change",
                   (yy ? "parameter" : "gradient"));
      return 3;
    }
    if (mp == m) --mp; /* oldest pair is lost */
    if (--mark < 0) mark = m - 1;
    ++isave(6); /* rejects */
    gamma = 1.0;
  }
  isave(9) = mp;  
  isave(10) = mark;
  //write, format=" (%d) MP=%d\n", mark, mp;
  if (mp < 1) return 0; /* Empty L-BFGS approximation: use steepest descent. */

  /********************************************************************
   ***                                                              ***
   ***  Obtain new search direction by the 2-loop L-BFGS recursion  ***
   ***                                                              ***
   ********************************************************************/

  alpha = array(double, m);
  k = mark + 1;
  if (constrained && (flags & OP_FLAG_RESTRICT_SOLVE)) {
    /* Restricted case: we must check every memorized pair. */
    gamma = 0.0;
    i = where(free); /* guaranteed to be non-empty */
    d = g(i);
    for (j=1 ; j<=mp ; ++j) {
      if (--k < 0) k = m - 1;
      if ((rho(k) = sum((s_k = (*s(k))(i))*(y_k = (*y(k))(i)))) > 0.0) {
        d -= (alpha(k) = sum(s_k*d)/rho(k))*y_k;
        if (! gamma) gamma = sum(y_k*y_k);
      }
      y_k = s_k = [];
    }
    if (gamma <= 0.0) return 0; /* use steepest descent */
    //if (OP_PRECOND == "auto" && mp >= 1) {
    //  if (1) {
    //    q0 = mark;
    //    q1 = (*s(q0))*(*y(q0));
    //    q2 = (*y(q0))*(*y(q0));
    //    for (j=2 ; j<=mp ; ++j) {
    //      if (--q0 < 0) q0 = m - 1;
    //      q1 += (*s(q0))*(*y(q0));
    //      q2 += (*y(q0))*(*y(q0));
    //    }
    //    q3 = array(gamma, dimsof(d));
    //    q3(i) = q1(i)/q2(i);
    // } else {
    //    q0 = mark;
    //    q1 = (*s(q0))*(*s(q0));
    //    q2 = (*y(q0))*(*y(q0));
    //    for (j=2 ; j<=mp ; ++j) {
    //      if (--q0 < 0) q0 = m - 1;
    //      q1 += (*s(q0))*(*s(q0));
    //      q2 += (*y(q0))*(*y(q0));
    //    }
    //    q3 = array(gamma, dimsof(d));
    //    q3(i) = sqrt(q1(i)/q2(i));
    //  }
    //  fma;pli,q3;pause,1;
    //  d = q3*d;
    //} else if (gamma != 1) {
    //  d = gamma*unref(d);
    //}
    if (gamma != 1) d = gamma*unref(d); /* FIXME: pre-conditioning */
    for (j=1 ; j<=mp ; ++j) {
      if (rho(k) > 0.0) {
        d += (alpha(k) - sum((*y(k))(i)*d)/rho(k))*(*s(k))(i);
      }
      if (++k >= m) k = 0;
    }
    (tmp = array(double, dimsof(x)))(i) = unref(d);
    eq_nocopy, d, tmp;
    i = []; /* save memory */
  } else {
    /* Unrestricted case: we already know that all memorized pairs are
       valid (i.e. rho(k) > zero). */
    eq_nocopy, d, g;
    for (j=1 ; j<=mp ; ++j) {
      if (--k < 0) k = m - 1;
      d -= (alpha(k) = sum((*s(k))*d)/rho(k))*(*y(k));
    }
    //if (OP_PRECOND == "auto" && mp >= 1) {
    //  if (1) {
    //    q0 = mark;
    //    q1 = (*s(q0))*(*y(q0));
    //    q2 = (*y(q0))*(*y(q0));
    //    for (j=2 ; j<=mp ; ++j) {
    //      if (--q0 < 0) q0 = m - 1;
    //      q1 += (*s(q0))*(*y(q0));
    //      q2 += (*y(q0))*(*y(q0));
    //    }
    //    q3 = q1/q2;
    //  } else {
    //    q0 = mark;
    //    q1 = (*s(q0))*(*s(q0));
    //    q2 = (*y(q0))*(*y(q0));
    //    for (j=2 ; j<=mp ; ++j) {
    //      if (--q0 < 0) q0 = m - 1;
    //      q1 += (*s(q0))*(*s(q0));
    //      q2 += (*y(q0))*(*y(q0));
    //    }
    //    q3 = sqrt(q1/q2);
    //  }
    //  fma;pli,q3;pause,1;
    //  d = q3*d;
    //} else if (gamma != 1) {
    //  d = gamma*unref(d);
    //}
    if (gamma != 1) d = gamma*unref(d); /* FIXME: pre-conditioning */
    for (j=1 ; j<=mp ; ++j) {
      d += (alpha(k) - sum((*y(k))*d)/rho(k))*(*s(k));
      if (++k >= m) k = 0;
    }
  }

  if ((dg = -sum(d*g)) >= 0.0) {
     /* Recursion failed to provide a descent direction: restart algorithm. */
    return 0;
  }
  return 1;
}

/*---------------------------------------------------------------------------*/
/* LINE SEARCH ROUTINES */

func op_csrch(f, g, &stp, sftol, sgtol, sxtol, stpmin, stpmax, &task,
              isave, dsave)
/* DOCUMENT op_csrch, f, g, stp, sftol, sgtol, sxtol, stpmin, stpmax, task,
 *                    isave, dsave;
 *
 *   This subroutine  finds a  step that  satisfies  a  sufficient decrease
 *   condition and a curvature condition.
 *
 *   Each call of the subroutine updates an interval with endpoints STX and
 *   STY.  The interval is initially chosen so that it contains a minimizer
 *   of the modified function
 *
 *         psi(stp) = f(stp) - f(0) - sftol*stp*f'(0).
 *
 *   If psi(stp) <= 0 and f'(stp)  >= 0 for some step, then the interval is
 *   chosen  so  that  it contains  a  minimizer of  F.   The algorithm  is
 *   designed  to find  a  step  that  satisfies  the  sufficient  decrease
 *   condition
 *
 *         f(stp) <= f(0) + sftol*stp*f'(0),
 *
 *   and the curvature condition
 *
 *         abs(f'(stp)) <= sgtol*abs(f'(0)).
 *
 *   If  SFTOL is  less than  SGTOL and  if, for  example, the  function is
 *   bounded  below, then  there  is  always a  step  which satisfies  both
 *   conditions.  If no  step can be found that  satisfies both conditions,
 *   then  the algorithm  stops  with a  warning.   In this  case stp  only
 *   satisfies the sufficient decrease condition.
 *
 *   A typical invocation of dcsrch has the following outline:
 *
 *   |  task = "start";
 *   |  stp = ...; // Initial STP value
 *   |  for (;;) {
 *   |    if (task == "fg" || task == "start") {
 *   |      // Evaluate the function and the gradient at STP.
 *   |    } else if (strpart(task, 1:4) == "conv") {
 *   |      // Search has converged.
 *   |    } else if (strpart(task, 1:4) == "warn") {
 *   |      // Some problem prevents further progress.
 *   |    } else {
 *   |      // An error occured.
 *   |    }
 *   |    op_csrch, ...;
 *   |  }
 *
 *   NOTE: The user must not alter work arrays between calls.
 *
 *   The subroutine arguments are:
 *
 *     F is a double precision  variable.  On initial entry, F is the value
 *       of  the function at 0.  On  subsequent entries, F is  the value of
 *       the function at STP.  On exit, F is left unchanged.
 *
 *     G is  a double  precision  variable.  On  initial entry,  G  is  the
 *       derivative of the function at  0.  On subsequent entries, G is the
 *       derivative of the function at STP.  On exit, G is left unchanged.
 *
 *    STP is a double  precision  variable.  On entry, STP  is the  current
 *       estimate  of a  satisfactory step.  On  initial entry,  a positive
 *       initial  estimate must be  provided.   On exit with TASK="fg", STP
 *       is  the  new  estimate  of a  satisfactory  step.   On  exit  with
 *       strpart(TASK,1:3)="conv",   STP is  left  unchanged and  satisfies
 *       the sufficient  decrease and curvature  condition.   On exit  with
 *       TASK not equal to "fg", STP is left unchanged.
 *
 *     SFTOL is a  double precision variable.  On entry,  SFTOL specifies a
 *       nonnegative tolerance  for the sufficient  decrease condition.  On
 *       exit, SFTOL is unchanged.
 *
 *     SGTOL is a  double precision variable.  On entry,  SGTOL specifies a
 *       nonnegative tolerance for the curvature condition.  On exit, SGTOL
 *       is unchanged.
 *
 *     SXTOL is a  double precision variable.  On entry,  SXTOL specifies a
 *       nonnegative  relative  tolerance  for  an  acceptable  step.   The
 *       subroutine exits with a warning if the relative difference between
 *       STY and STX is less than SXTOL.  On exit, SXTOL is unchanged.
 *
 *     STPMIN  is  a  double precision  variable.   On entry,  STPMIN is  a
 *       nonnegative  lower  bound  for  the step.   On  exit,  STPMIN   is
 *       unchanged.
 *
 *     STPMAX is  a  double precision   variable.   On entry,  STPMAX is  a
 *       nonnegative   upper  bound  for  the step.   On  exit,  STPMAX  is
 *       unchanged.
 *
 *     TASK is a character variable of length at least 60.
 *       On initial entry, task must be set to "start".
 *       On exit, TASK indicates the required action:
 *
 *          If task(1:2) = "FG" then evaluate the function and
 *          derivative at stp and call dcsrch again.
 *
 *          If task(1:4) = "CONV" then the search is successful.
 *
 *          If task(1:4) = "WARN" then the subroutine is not able
 *          to satisfy the convergence conditions. The exit value of
 *          stp contains the best point found during the search.
 *
 *          If task(1:5) = "ERROR" then there is an error in the
 *          input arguments.
 *
 *       On exit with convergence, a warning or an error, the variable TASK
 *       contains additional information.
 *
 *     ISAVE is an integer work array of, at least, 2 elements.
 *
 *     DSAVE is a double precision work array of, at least, 13 elements.
 *
 *
 * HISTORY:
 *   MINPACK-1 Project. June 1983.
 *   Argonne National Laboratory.
 *   Jorge J. More' and David J. Thuente.
 *
 *   MINPACK-2 Project. November 1993.
 *   Argonne National Laboratory and University of Minnesota.
 *   Brett M. Averick, Richard G. Carter, and Jorge J. More'.
 *
 *   Yorick translation an improvements.  October 2001.
 *   Observatoire de Lyon.
 *   Eric Thiebaut.
 *
 *
 * SEE ALSO: op_cstep, op_cgmn, op_vmlm.
 */
{
  /* Initialization block.*/
  zero = 0.0;
  xtrapl = 1.1;
  xtrapu = 4.0;

  //write, format="%s: SF_TOL=%g SG_TOL=%g SX_TOL=%g\n", task, sftol, sgtol, sxtol;
  if (strpart(task, 1:5) == "start") {

    /* Check the input arguments for errors.
       Exit if there are errors on input. */
    if (stpmax < stpmin) { task = "error: stpmax < stpmin";    return; }
    if (stpmin < zero)   { task = "error: stpmin < zero";      return; }
    if (sxtol < zero)    { task = "error: sxtol < zero";       return; }
    if (sgtol < zero)    { task = "error: sgtol < zero";       return; }
    if (sftol < zero)    { task = "error: sftol < zero";       return; }
    if (g >= zero)       { task = "error: initial g >= zero";  return; }
    if (stp > stpmax)    { task = "error: stp > stpmax";       return; }
    if (stp < stpmin)    { task = "error: stp < stpmin";       return; }

    /* Initialize local variables.
       The variables stx, fx, gx contain the values of the step,
       function, and derivative at the best step.
       The variables sty, fy, gy contain the value of the step,
       function, and derivative at sty.
       The variables stp, f, g contain the values of the step,
       function, and derivative at stp. */
    brackt = 0;
    stage = 1;
    finit = f;
    ginit = g;
    gtest = sftol*ginit;
    width = stpmax - stpmin;
    width1 = 2.0*width;
    stx = zero;
    fx = finit;
    gx = ginit;
    sty = zero;
    fy = finit;
    gy = ginit;
    stmin = zero;
    stmax = stp + xtrapu*stp;
    task = "fg";

  } else {

    /* Restore local variables. */
    brackt = isave( 1);
    stage  = isave( 2);
    ginit  = dsave( 1);
    gtest  = dsave( 2);
    gx     = dsave( 3);
    gy     = dsave( 4);
    finit  = dsave( 5);
    fx     = dsave( 6);
    fy     = dsave( 7);
    stx    = dsave( 8);
    sty    = dsave( 9);
    stmin  = dsave(10);
    stmax  = dsave(11);
    width  = dsave(12);
    width1 = dsave(13);

    /* If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the
       algorithm enters the second stage.*/
    ftest = finit + stp*gtest;
    if (stage == 1 && f <= ftest && g >= zero) stage = 2;

    /* Test for termination: convergence or warnings. */
    if (f <= ftest && abs(g) <= sgtol*(-ginit)) {
      task = "convergence";
    } else if (stp == stpmin && (f > ftest || g >= gtest)) {
      task = "warning: stp = stpmin";
    } else if (stp == stpmax && f <= ftest && g <= gtest) {
      task = "warning: stp = stpmax";
    } else if (brackt && stmax - stmin <= sxtol*stmax) {
      task = "warning: sxtol test satisfied";
    } else if (brackt && (stp <= stmin || stp >= stmax)) {
      task = "warning: rounding errors prevent progress";
    } else {
      /* A modified function is used to predict the step during the first
         stage if a lower function value has been obtained but the decrease
         is not sufficient. */
      if (stage == 1 && f <= fx && f > ftest) {
        /* Define the modified function and derivative values. */
        fm  = f  - stp*gtest;
        fxm = fx - stx*gtest;
        fym = fy - sty*gtest;
        gm  = g  - gtest;
        gxm = gx - gtest;
        gym = gy - gtest;

        /* Call dcstep to update stx, sty, and to compute the new step. */
        op_cstep,stx,fxm,gxm,sty,fym,gym,stp,fm,gm,brackt,stmin,stmax;

        /* Reset the function and derivative values for f. */
        fx = fxm + stx*gtest;
        fy = fym + sty*gtest;
        gx = gxm + gtest;
        gy = gym + gtest;

      } else {
        /* Call dcstep to update stx, sty, and to compute the new step. */
        op_cstep,stx,fx,gx,sty,fy,gy,stp,f,g,brackt,stmin,stmax;
      }

      if (brackt) {
        /* Decide if a bisection step is needed. */
        if ((w = abs(sty - stx)) >= 0.66*width1) stp = stx + 0.5*(sty - stx);
        width1 = width;
        width = w;

        /* Set the minimum and maximum steps allowed for stp. */
        stmin = min(stx, sty);
        stmax = max(stx, sty);
      } else {
        /* Set the minimum and maximum steps allowed for stp. */
        stmin = stp + xtrapl*(stp - stx);
        stmax = stp + xtrapu*(stp - stx);
      }

      /* Force the step to be within the bounds stpmax and stpmin. */
      stp = min(max(stp, stpmin), stpmax);

      /* If further progress is not possible, let stp be the best
         point obtained during the search. */
      if (brackt && (stp <= stmin || stp >= stmax)
          || (brackt && stmax-stmin <= sxtol*stmax)) stp = stx;

      /* Obtain another function and derivative. */
      task = "fg";
    }
  }

  /* Save local variables. */
  isave( 1) = brackt;
  isave( 2) = stage;
  dsave( 1) = ginit;
  dsave( 2) = gtest;
  dsave( 3) = gx;
  dsave( 4) = gy;
  dsave( 5) = finit;
  dsave( 6) = fx;
  dsave( 7) = fy;
  dsave( 8) = stx;
  dsave( 9) = sty;
  dsave(10) = stmin;
  dsave(11) = stmax;
  dsave(12) = width;
  dsave(13) = width1;
}

func op_cstep(&stx, &fx, &dx, &sty, &fy, &dy, &stp, fp, dp,
              &brackt, stpmin, stpmax)
/* DOCUMENT op_cstep, stx, fx, dx,
 *                    sty, fy, dy,
 *                    stp, fp, dp,
 *                    brackt, stpmin, stpmax;
 *
 *   This subroutine computes a safeguarded step for a search procedure and
 *   updates an  interval that contains a step  that satisfies a sufficient
 *   decrease and a curvature condition.
 *
 *   The parameter STX contains the  step with the least function value. If
 *   brackt  is set  to .true. then  a minimizer  has been bracketed  in an
 *   interval with  endpoints STX and sty.  The  parameter STP contains the
 *   current step.  The subroutine assumes that if BRACKT is true then:
 *
 *         min(STX,STY) < STP < max(STX,STY),
 *
 *   and  that the derivative  at STX is  negative in the direction  of the
 *   step.
 *
 *   The subroutine arguments are:
 *
 *     STX is a double precision  variable.  On entry, STX is the best step
 *       obtained so  far and is an endpoint of  the interval that contains
 *       the minimizer.  On exit, STX is the updated best step.
 *
 *     FX is a double precision variable.  On entry, FX  is the function at
 *       STX.  On exit, FX is the function at STX.
 *
 *     DX is  a double precision variable.  On entry,  DX is the derivative
 *       of  the function at STX.   The derivative must be  negative in the
 *       direction  of  the  step,  that is,  DX  and  STP - STX must  have
 *       opposite signs.  On exit, DX  is the derivative of the function at
 *       STX.
 *
 *     STY  is a double  precision variable.  On  entry, STY is  the second
 *       endpoint  of the interval  that contains the minimizer.   On exit,
 *       STY is the  updated  endpoint  of the  interval that  contains the
 *       minimizer.
 *
 *     FY is a double precision  variable.  On entry, FY is the function at
 *       STY.  On exit, FY is the function at STY.
 *
 *     DY is a double precision variable.
 *       On entry, DY is the derivative of the function at STY.
 *       On exit, DY is the derivative of the function at the exit STY.
 *
 *     STP  is a double precision  variable.  On entry, STP  is the current
 *       step. If  BRACKT is true, then, on input, STP  must be between STX
 *       and STY.  On exit, STP is a new trial step.
 *
 *     FP is a double precision  variable.  On entry, FP is the function at
 *       STP.  On exit, FP is unchanged.
 *
 *     DP is  a double precision variable.  On entry,  DP is the derivative
 *       of the function at STP.  On exit, DP is unchanged.
 *
 *     BRACKT is  a logical  variable.  On  entry, BRACKT  specifies  if  a
 *       minimizer  has been  bracketed.  Initially BRACKT  must be  set to
 *       FALSE.    On  exit,  BRACKT  specifies if  a  minimizer  has  been
 *       bracketed.  When a minimizer is bracketed brackt is set to TRUE.
 *
 *     STPMIN is a double precision  variable.  On entry, STPMIN is a lower
 *       bound for the step.  On exit, STPMIN is unchanged.
 *
 *     STPMAX is a double precision variable.  On entry, STPMAX is an upper
 *       bound for the step.  On exit, STPMAX is unchanged.
 *
 *
 * HISTORY:
 *   MINPACK-1 Project. June 1983
 *   Argonne National Laboratory.
 *   Jorge J. More' and David J. Thuente.
 *
 *   MINPACK-2 Project. November 1993.
 *   Argonne National Laboratory and University of Minnesota.
 *   Brett M. Averick and Jorge J. More'.
 *
 *   Yorick translation an improvements.  October 2001.
 *   Observatoire de Lyon.
 *   Eric Thiebaut.
 *
 *
 * SEE ALSO: op_csrch.
 */
{
#if 0
  /* Check the input parameters for errors. */
  if ((brackt && (stx < sty ? (stp <= stx || stp >= sty)
		            : (stp >= stx || stp <= sty)))) {
    error, "STP outside bracket (STX,STY)";
  } else if (dx*(stp - stx) >= 0.0) {
    error, "descent condition violated";
  } else if (stpmax < stpmin) {
    error, "STPMAX < STPMIN";
  }
#endif

  /* Determine if the derivatives have opposite sign. */
  sgnd = dp*(dx/abs(dx));

  if (fp > fx) {
    /* First case: A higher function value. The minimum is bracketed.  If
       the cubic step is closer to stx than the quadratic step, the cubic
       step is taken, otherwise the average of the cubic and quadratic
       steps is taken. */
    theta = 3.0*(fx - fp)/(stp - stx) + dx + dp;
    s = max(abs(theta), abs(dx), abs(dp));
    gamma = s*sqrt((theta/s)^2 - (dx/s)*(dp/s));
    if (stp < stx) gamma = -gamma;
    p = (gamma - dx) + theta;
    q = ((gamma - dx) + gamma) + dp;
    r = p/q;
    stpc = stx + r*(stp - stx);
    stpq = stx + ((dx/((fx - fp)/(stp - stx) + dx))/2.0)*(stp - stx);
    if (abs(stpc - stx) < abs(stpq - stx)) {
      stpf = stpc;
    } else {
      /* stpf = (stpq + stpc)/2.0; */
      stpf = stpc + (stpq - stpc)/2.0;
    }
    brackt = 1;
  } else if (sgnd < 0.0) {
    /* Second case: A lower function value and derivatives of opposite
       sign. The minimum is bracketed. If the cubic step is farther from
       stp than the secant step, the cubic step is taken, otherwise the
       secant step is taken. */
    theta = 3.0*(fx - fp)/(stp - stx) + dx + dp;
    s = max(abs(theta),abs(dx),abs(dp));
    gamma = s*sqrt((theta/s)^2 - (dx/s)*(dp/s));
    if (stp > stx) gamma = -gamma;
    p = (gamma - dp) + theta;
    q = ((gamma - dp) + gamma) + dx;
    r = p/q;
    stpc = stp + r*(stx - stp);
    stpq = stp + (dp/(dp - dx))*(stx - stp);
    if (abs(stpc - stp) > abs(stpq - stp)) {
      stpf = stpc;
    } else {
      stpf = stpq;
    }
    brackt = 1;
  } else if (abs(dp) < abs(dx)) {
    /* Third case: A lower function value, derivatives of the same sign,
       and the magnitude of the derivative decreases.  The cubic step is
       computed only if the cubic tends to infinity in the direction of the
       step or if the minimum of the cubic is beyond stp. Otherwise the
       cubic step is defined to be the secant step. */
    theta = 3.0*(fx - fp)/(stp - stx) + dx + dp;
    s = max(abs(theta),abs(dx),abs(dp));
    /* The case gamma = 0 only arises if the cubic does not tend to
       infinity in the direction of the step. */
    gamma = s*sqrt(max(0.0, (theta/s)^2 - (dx/s)*(dp/s)));
    if (stp > stx) gamma = -gamma;
    p = (gamma - dp) + theta;
    q = (gamma + (dx - dp)) + gamma;
    r = p/q;
    if (r < 0.0 && gamma != 0.0) {
      stpc = stp + r*(stx - stp);
    } else if (stp > stx) {
      stpc = stpmax;
    } else {
      stpc = stpmin;
    }
    stpq = stp + (dp/(dp - dx))*(stx - stp);
    if (brackt) {
      /* A minimizer has been bracketed. If the cubic step is closer to stp
	 than the secant step, the cubic step is taken, otherwise the
	 secant step is taken. */
      if (abs(stpc - stp) < abs(stpq - stp)) {
	stpf = stpc;
      } else {
	stpf = stpq;
      }
      if (stp > stx) {
	stpf = min(stp + 0.66*(sty - stp), stpf);
      } else {
	stpf = max(stp + 0.66*(sty - stp), stpf);
      }
    } else {
      /* A minimizer has not been bracketed. If the cubic step is farther
	 from stp than the secant step, the cubic step is taken, otherwise
	 the secant step is taken. */
      if (abs(stpc-stp) > abs(stpq-stp)) {
	stpf = stpc;
      } else {
	stpf = stpq;
      }
      stpf = min(stpmax,stpf);
      stpf = max(stpmin,stpf);
    }
  } else {
    /* Fourth case: A lower function value, derivatives of the same sign,
       and the magnitude of the derivative does not decrease. If the
       minimum is not bracketed, the step is either stpmin or stpmax,
       otherwise the cubic step is taken. */
    if (brackt) {
      theta = 3.0*(fp - fy)/(sty - stp) + dy + dp;
      s = max(abs(theta),abs(dy),abs(dp));
      gamma = s*sqrt((theta/s)^2 - (dy/s)*(dp/s));
      if (stp > sty) gamma = -gamma;
      p = (gamma - dp) + theta;
      q = ((gamma - dp) + gamma) + dy;
      r = p/q;
      stpc = stp + r*(sty - stp);
      stpf = stpc;
    } else if (stp > stx) {
      stpf = stpmax;
    } else {
      stpf = stpmin;
    }
  }

  /* Update the interval which contains a minimizer. */
  if (fp > fx) {
    sty = stp;
    fy = fp;
    dy = dp;
  } else {
    if (sgnd < 0.0) {
      sty = stx;
      fy = fx;
      dy = dx;
    }
    stx = stp;
    fx = fp;
    dx = dp;
  }

  /* Compute the new step. */
  stp = stpf;
}

/*---------------------------------------------------------------------------*/
/* UTILITIES */

func op_not_scalar_integer(&val, def) {
  if (is_void(val)) val = def;
  if (! is_array(val) || dimsof(val)(1)) return 1; /* not a scalar */
  if ((s = structof(val)) == long) return 0;
  if (s==int || s==short || s==char) { val = long(val); return 0; }
  return 1; /* bada data type */
}
func op_not_scalar_real(&val, def)
/* DOCUMENT op_not_scalar_integer(var, def)
       -or- op_not_scalar_real(var, def)

     Make sure  contents of symbol VAR  is a scalar integer  (long) or real
     (double).  If  VAR is nil,  then it is  set to the default  value DEF.
     VAR is converted to long/double as needed.  If result is true then VAR
     (or DEF) is  not a valid scalar integer/real.   For instance, inside a
     function where  OPT is  an optional keyword/parameter  that must  be a
     scalar real:

        if (op_not_scalar_real(opt, 0.0)) error, "bad value for OPT";

   SEE ALSO: */
{
  if (is_void(val)) val = def;
  if (! is_array(val) || dimsof(val)(1)) return 1; /* not a scalar */
  if ((s = structof(val)) == double) return 0;
  if (s==float || s==long || s==int || s==short || s==char) {
    val = double(val);
    return 0;
  }
  return 1;
}

func op_unref(&x)
/* DOCUMENT op_unref(x)
     returns X, destroying X in the process (useful to deal with temporary
     big arrays).  Written after Yorick's FAQ.
   SEE ALSO: eq_nocopy, unref, swap. */
{
  local y;
  eq_nocopy, y, x;
  x = [];
  return y;
}

func op_nrm2(x)
/* DOCUMENT op_nrm2(x)
    Returns the Euclidian norm of X: sqrt(sum(X*X)), taking care of
    overflows.  X can be complex.
*/
{
  if (structof(x) == complex) {
    y = x.im;
    x = double(x);
    if (! (s = max(-min(x), max(x), -min(y), max(y)))) return 0.0;
    if (s == 1.0) return sqrt(sum(x*x) + sum(y*y));
    q = 1.0/s;
    x *= q;
    y *= q;
    return s*sqrt(sum(x*x) + sum(y*y));
  } else {
    if (! (s = max(-min(x), max(x)))) return 0.0;
    if (s == 1.0) return sqrt(sum(x*x));
    x *= 1.0/s;
    return s*sqrt(sum(x*x));
  }
}

func op_apply_bounds(xmin, x, xmax)
/* DOCUMENT op_apply_bounds(xmin, x, xmax)
     Returns X after applying lower and/or upper bound constraints.  XMIN
     and/or XMAX can be nil if there is no such bound.
   SEE ALSO op_project_gradient. */
{
  if (is_void(xmin)) {
    if (is_void(xmax)) return x;
    return min(x, xmax);
  }
  if (is_void(xmax)) return max(x, xmin);
  return min(max(x, xmin), xmax);
}

func op_project_gradient(g, xmin, x, xmax)
/* DOCUMENT op_project_gradient(g, xmin, x, xmax)
     Returns projected gradient G according to bounds XMIN and/or XMAX for
     X.  XMIN and/or XMAX can be nil if there is no such bound.
   SEE ALSO op_apply_bounds. */
{
  if (is_void(xmin)) {
    if (is_void(xmax)) return g;
    return ((x >= xmax) & (g < 0.0))*g;
  }
  if (is_void(xmax)) return ((x <= xmin) & (g > 0.0))*g;
  return (((x >= xmax) & (g < 0.0)) | ((x <= xmin) & (g > 0.0)))*g;
}

/*
 * Forward search, constraints read:
 *           XMIN <= X0 + STP*D <= XMAX
 *     <==>  XMIN - X0 <= STP*D <= XMAX - X0
 *     <==>  STP <= (XMAX - X0)/D     where D > 0
 *           STP <= (XMIN - X0)/D     where D < 0
 *
 * Backward search, constraints read:
 *           XMIN <= X0 - STP*D <= XMAX
 *     <==>  X0 - XMAX <= STP*D <= X0 - XMIN
 *     <==>  STP <= (X0 - XMIN)/D     where D > 0
 *           STP <= (X0 - XMAX)/D     where D < 0
 */
func op_get_max_step(stpmax, d, xmin, x0, xmax)
{
  sbig = abs(stpmax);
  found = 0n;
  if (stpmax > 0.0) {
    /* Forward search. */
    eq_nocopy, d, double(d); /* a no-op if D is already of that type */
    if (! is_void(xmin) && is_array((i = where(d < 0.0)))) {
      smax = max((xmin - x0)(i)/d(i));
      found = 1n;
    } else {
      smax = -sbig;
    }
    if (! is_void(xmax) && is_array((i = where(d > 0.0)))) {
      smax = max(smax, max((xmax - x0)(i)/d(i)));
      found = 1n;
    }
  } else {
    /* Backward search. */
    eq_nocopy, d, double(d); /* a no-op if D is already of that type */
    if (! is_void(xmin) && is_array((i = where(d > 0.0)))) {
      smax = max((x0 - xmin)(i)/d(i));
      found = 1n;
    } else {
      smax = -sbig;
    }
    if (! is_void(xmax) && is_array((i = where(d < 0.0)))) {
      smax = max(smax, max((x0 - xmax)(i)/d(i)));
      found = 1n;
    }
  }
  return (found ? smax : sbig);
}

func op_get_step_bounds(stpmax, d, xmin, x0, xmax)
/* DOCUMENT [stp1, stp2] = op_get_step_bounds(stpmax, d, xmin, x0, xmax)
      Compute  particular step  values  for bound  constrained line  search
      along direction  D.  If STPMAX >  0, then D is  the search direction;
      otherwise,  then  D  is  the  anti-search direction.   In  any  case,
      abs(STPMAX) is  the maximum step  size to consider.  In  other words,
      starting  from X0,  the  line search  seeks  for an  improved set  of
      parameters X such that:

         X = X0 + STP*D     if STPMAX > 0
         X = X0 - STP*D     if STPMAX < 0

      where the _positive_ step STP is in the range [0, abs(STPMAX)].

      On return, STP1 is the step to the closest bound (i.e.  no bounds are
      encountered along  D for  steps smaller than  STP1), and STP2  is the
      step to  the farthest bound (i.e.  it is pointless to  seek for steps
      greater  than STP2).   STP1  and  STP2 are  in  [0, abs(STPMAX)]  and
      STP1=STP2=abs(STPMAX)  if no bounds  are encountered  along D  in the
      allowed range.  If, on return, STP1 < 0, then either input X violates
      some of the bounds or bounds  are incompatible (i.e. a lower bound is
      greater than an upper bound).

      XMIN and XMIN are the bounds on X  and may be nil if there is no such
      bound.
 */
{
  /* Compute step bounds. */
  smin = (sbig = abs(stpmax));
  smax = -sbig;
  if (stpmax > 0.0) {
    /* Forward search. */
    eq_nocopy, d, double(d); /* a no-op if D is already of that type */
    if (! is_void(xmin) && is_array((i = where(d < 0.0)))) {
      s = (xmin - x0)(i)/d(i);
      smin = min(smin, min(s));
      smax = max(smax, max(s));
      s = [];
    }
    if (! is_void(xmax) && is_array((i = where(d > 0.0)))) {
      s = (xmax - x0)(i)/d(i);
      smin = min(smin, min(s));
      smax = max(smax, max(s));
      s = [];
    }
  } else {
    /* Backward search. */
    eq_nocopy, d, double(d); /* a no-op if D is already of that type */
    if (! is_void(xmin) && is_array((i = where(d > 0.0)))) {
      s = (x0 - xmin)(i)/d(i);
      smin = min(smin, min(s));
      smax = max(smax, max(s));
      s = [];
    }
    if (! is_void(xmax) && is_array((i = where(d < 0.0)))) {
      s = (x0 - xmax)(i)/d(i);
      smin = min(smin, min(s));
      smax = max(smax, max(s));
      s = [];
    }
  }

  /* Return bounds (take care of the unbounded case for SMAX). */
  return [smin, (smax >= smin ? smax : sbig)];
}

/*---------------------------------------------------------------------------*/
/* DEALING WITH INTERPRETED/BUILT-IN VERSION OF OPTIMPACK ROUTINES */

func _op_restore_builtin(builtin, interpreted, name)
{
  if (is_func(builtin) == 2) return builtin;
  if (is_func(interpreted) == 1){
    write, format="warning: only got interpreted version of \"%s\"\n", name;
    return interpreted;
  }
  write, format="warning: no function \"%s\"\n", name; 
}

func _op_restore_interpreted(builtin, interpreted, name)
{
  if (is_func(interpreted) == 1) return interpreted;
  if (is_func(builtin) == 2) {
    write, format="warning: only got built-in version of \"%s\"\n", name;
    return builtin;
  }
  write, format="warning: no function \"%s\"\n", name; 
}

local OP_RESTORE_CODE;
func op_restore(what)
/* DOCUMENT op_restore, what;
     Restores interpreted or built-in version of (some) OptimPack functions.
     Argument WHAT is a string:
       "auto"        - restore built-in code if any, otherwise interpreted;
       "builtin"     - restore built-in version;
       "interpreted" - restore interpreted version.
     This function is needed to help me having only one Yorick file with
     test code and documentation of OptimPack interpreted/built-in
     routines.  The default (when optimpack.i is parsed the first time by
     Yorick) is to attempt to use built-in version of the functions
     (WHAT="auto").
 */
{
  if (what == "auto") {
    restore = (is_func(op_vmlmb_next_BUILTIN) == 2 ?
               _op_restore_builtin :
               _op_restore_interpreted);
  } else if (what == "builtin") {
    restore = _op_restore_builtin;
  } else if (what == "interpreted") {
    restore = _op_restore_interpreted;
  } else {
    error, "syntax: op_setup, \"auto|builtin|interpreted\"";
  }
  {extern OP_RESTORE_CODE; OP_RESTORE_CODE = what;}

  /* Restore any saved function (these statements are produced by the
     script optimpack-fix). */
  {extern op_cgmnb_next;}
  op_cgmnb_next = restore(op_cgmnb_next_BUILTIN, op_cgmnb_next_INTERPR,
                          "op_cgmnb_next");
  {extern op_cgmnb_setup;}
  op_cgmnb_setup = restore(op_cgmnb_setup_BUILTIN, op_cgmnb_setup_INTERPR,
                           "op_cgmnb_setup");
  {extern op_get_amaxg;}
  op_get_amaxg = restore(op_get_amaxg_BUILTIN, op_get_amaxg_INTERPR,
                         "op_get_amaxg");
  {extern op_get_flags;}
  op_get_flags = restore(op_get_flags_BUILTIN, op_get_flags_INTERPR,
                         "op_get_flags");
  {extern op_get_iter;}
  op_get_iter = restore(op_get_iter_BUILTIN, op_get_iter_INTERPR,
                        "op_get_iter");
  {extern op_get_msg;}
  op_get_msg = restore(op_get_msg_BUILTIN, op_get_msg_INTERPR,
                       "op_get_msg");
  {extern op_get_rejects;}
  op_get_rejects = restore(op_get_rejects_BUILTIN, op_get_rejects_INTERPR,
                           "op_get_rejects");
  {extern op_get_restarts;}
  op_get_restarts = restore(op_get_restarts_BUILTIN, op_get_restarts_INTERPR,
                            "op_get_restarts");
  {extern op_get_stage;}
  op_get_stage = restore(op_get_stage_BUILTIN, op_get_stage_INTERPR,
                         "op_get_stage");
  {extern op_get_step;}
  op_get_step = restore(op_get_step_BUILTIN, op_get_step_INTERPR,
                        "op_get_step");
  {extern op_vmlmb_next;}
  op_vmlmb_next = restore(op_vmlmb_next_BUILTIN, op_vmlmb_next_INTERPR,
                          "op_vmlmb_next");
  {extern op_vmlmb_setup;}
  op_vmlmb_setup = restore(op_vmlmb_setup_BUILTIN, op_vmlmb_setup_INTERPR,
                           "op_vmlmb_setup");
}

/*---------------------------------------------------------------------------*/
/* Save interpreted version of functions implemented by OptimPack
   (these statements _must_ be at the end of this file and are
   produced by the script optimpack-fix); then restore either built-in
   or interpreted functions. */

op_cgmnb_next_INTERPR = op_cgmnb_next;
op_cgmnb_setup_INTERPR = op_cgmnb_setup;
op_get_amaxg_INTERPR = op_get_amaxg;
op_get_flags_INTERPR = op_get_flags;
op_get_iter_INTERPR = op_get_iter;
op_get_msg_INTERPR = op_get_msg;
op_get_rejects_INTERPR = op_get_rejects;
op_get_restarts_INTERPR = op_get_restarts;
op_get_stage_INTERPR = op_get_stage;
op_get_step_INTERPR = op_get_step;
op_vmlmb_next_INTERPR = op_vmlmb_next;
op_vmlmb_setup_INTERPR = op_vmlmb_setup;

if (structof(OP_RESTORE_CODE) != string || dimsof(OP_RESTORE_CODE)(1) ||
    (OP_RESTORE_CODE != "auto" && OP_RESTORE_CODE != "builtin" &&
     OP_RESTORE_CODE != "interpreted")) OP_RESTORE_CODE = "auto";
op_restore, OP_RESTORE_CODE;

/*---------------------------------------------------------------------------*/
