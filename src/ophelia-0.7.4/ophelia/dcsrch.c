/**
 * dcsrch.c
 *
 * NLP (NonLinear Programming) Library
 * Copyright (C) 2002-2006 Laurent Martin
 *                         Emilie Wanufelle <ewa@math.fundp.ac.be>,
 *                         Fabian Bastin    <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

/**
 * The following code is adapted from the original routines of More
 * and Thuente, Minpack Project, June 1983.
 * We also refer to L-BFGS-B of Zhu, Byrd and Nocedal, 1994.
 *
 * The adaptation to C has been done by Laurent Martin and Emilie
 * Wanufelle, under the supervision of Fabian Bastin.
 *
 * Moreover, we also offer the possibility to compute both function and
 * gradient with the same function.
 */

#include "linesearch.h"

#include <stdio.h>

#include <oratio/ntools.h>

int dcstep(double  *stx, double *fx, double *dx, double *sty, double *fy,  
	   double *dy, double *stp, double fp, double dp, int *brackt, 
	   double stpmin, double stpmax, int *info)
{
  /**
   *     subroutine dcstep 
   
   *     the purpose of dcstep is to compute a safeguarded step for 
   *     a linesearch and to update an interval of uncertainty for 
   *     a minimizer of the function. 

   *     the parameter stx contains the step with the least function 
   *     value. the parameter stp contains the current step. it is 
   *     assumed that the derivative at stx is negative in the 
   *     direction of the step. if brackt is set true then a 
   *     minimizer has been bracketed in an interval of uncertainty 
   *     with endpoints stx and sty. 

   *     the subroutine statement is 

   *       subroutine dcstep(stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt, 
   *                        stpmin,stpmax,info) 

   *     where 

   *       stx, fx, and dx are variables which specify the step, 
   *         the function, and the derivative at the best step obtained 
   *         so far. the derivative must be negative in the direction 
   *         of the step, that is, dx and stp-stx must have opposite 
   *         signs. on output these parameters are updated appropriately. 

   *       sty, fy, and dy are variables which specify the step, 
   *         the function, and the derivative at the other endpoint of 
   *         the interval of uncertainty. on output these parameters are 
   *         updated appropriately. 
   
   *       stp, fp, and dp are variables which specify the step, 
   *         the function, and the derivative at the current step. 
   *         if brackt is set true then on input stp must be 
   *         between stx and sty. on output stp is set to the new step. 

   *       brackt is a long int variable which specifies if a minimizer 
   *         has been bracketed. if the minimizer has not been bracketed 
   *         then on input brackt must be set false. if the minimizer 
   *         is bracketed then on output brackt is set true. 

   *       stpmin and stpmax are input variables which specify lower 
   *         and upper bounds for the step. 

   *       info is an int output variable set as follows: 
   *         if info = 1,2,3,4,5, then the step has been computed 
   *         according to one of the five cases below. otherwise 
   *         info = 0, and this indicates improper input parameters. 

   *     subprograms called 

   *     argonne national laboratory. minpack project. june 1983 
   *     jorge j. more', david j. thuente 
   */

  /* system generated locals */
  double d_1, d_2, d_3;

  /* local variables */
  static double sgnd, stpc, stpf, stpq, p, q, gamma, r, s, theta;
  static int bound;

  *info = 0;

  /* check the input parameters for errors. */
  if ( ( *brackt && 
	 (*stp <= nt_min(*stx,*sty) || *stp >= nt_max(*stx,*sty)) ) || 
       *dx * (*stp - *stx) >= 0. || stpmax < stpmin) {
    return 0;
  }

  /* determine if the derivatives have opposite sign. */
    
  sgnd = dp * (*dx / fabs(*dx));
    
  /* first case. a higher function value. */
  /* the minimum is bracketed. if the cubic step is closer */
  /* to stx than the quadratic step, the cubic step is taken, */
  /* else the average of the cubic and quadratic steps is taken. */
    
  if (fp > *fx || isnan(fp)) {
    *info = 1;
    bound = 1;
    theta = (*fx - fp) * 3 / (*stp - *stx) + *dx + dp;
    /* computing max */
    d_1 = fabs(theta), d_2 = fabs(*dx), d_1 = nt_max(d_2,d_1), d_2 = fabs(dp);
    s = nt_max(d_2,d_1);
    /* computing 2nd power */
    d_1 = theta / s;
    gamma = s * sqrt(d_1 * d_1 - *dx / s * (dp / s));
    if (*stp < *stx) {
      gamma = -gamma;
    }
    p = gamma - *dx + theta;
    q = gamma - *dx + gamma + dp;
    r = p / q;
    stpc = *stx + r * (*stp - *stx);
    stpq = *stx + *dx / ((*fx - fp) / (*stp - *stx) + *dx) / 2 * (*stp - 
								  *stx);
    if ((d_1 = stpc - *stx, fabs(d_1)) < (d_2 = stpq - *stx, fabs(d_2))) {
      stpf = stpc;
    } else {
      stpf = stpc + (stpq - stpc) / 2;
    }
    *brackt = 1;
      
    /* second case. a lower function value and derivatives of */
    /* opposite sign. the minimum is bracketed. if the cubic */
    /* step is closer to stx than the quadratic (secant) step, */
    /* the cubic step is taken, else the quadratic step is taken. */
      
  } else if (sgnd < 0.) {
    *info = 2;
    bound = 0;
    theta = (*fx - fp) * 3 / (*stp - *stx) + *dx + dp;
    /* computing max */
    d_1 = fabs(theta), d_2 = fabs(*dx), d_1 = nt_max(d_2,d_1), d_2 = fabs(dp);
    s = nt_max(d_2,d_1);
    /* computing 2nd power */
    d_1 = theta / s;
    gamma = s * sqrt(d_1 * d_1 - *dx / s * (dp / s));
    if (*stp > *stx) {
      gamma = -gamma;
    }
    p = gamma - dp + theta;
    q = gamma - dp + gamma + *dx;
    r = p / q;
    stpc = *stp + r * (*stx - *stp);
    stpq = *stp + dp / (dp - *dx) * (*stx - *stp);
    if ((d_1 = stpc - *stp, fabs(d_1)) > (d_2 = stpq - *stp, fabs(d_2))) {
      stpf = stpc;
    } else {
      stpf = stpq;
    }
    *brackt = 1;
	
    /* third case. a lower function value, derivatives of the */
    /* same sign, and the magnitude of the derivative decreases. */
    /* the cubic step is only used if the cubic tends to infinity */
    /* in the direction of the step or if the minimum of the cubic */
    /* is beyond stp. otherwise the cubic step is defined to be */
    /* either stpmin or stpmax. the quadratic (secant) step is also */
	
    /* computed and if the minimum is bracketed then the the step */
    /* closest to stx is taken, else the step farthest away is taken. 
     */
	
  } else if (fabs(dp) < fabs(*dx)) {
    *info = 3;
    bound = 1;
    theta = (*fx - fp) * 3 / (*stp - *stx) + *dx + dp;
    /* computing max */
    d_1 = fabs(theta), d_2 = fabs(*dx), d_1 = nt_max(d_2,d_1), d_2 = fabs(dp);
    s = nt_max(d_2,d_1);
	
    /* the case gamma = 0 only arises if the cubic does not tend */
    /* to infinity in the direction of the step. */
	
    /* computing max */
    /* computing 2nd power */
    d_3 = theta / s;
    d_1 = 0., d_2 = d_3 * d_3 - *dx / s * (dp / s);
    gamma = s * sqrt((nt_max(d_2,d_1)));
    if (*stp > *stx) {
      gamma = -gamma;
    }
    p = gamma - dp + theta;
    q = gamma + (*dx - dp) + gamma;
    r = p / q;
    if (r < 0. && gamma != 0.) {
      stpc = *stp + r * (*stx - *stp);
    } else if (*stp > *stx) {
      stpc = stpmax;
    } else {
      stpc = stpmin;
    }
    stpq = *stp + dp / (dp - *dx) * (*stx - *stp);
    if (*brackt) {
      if ((d_1 = *stp - stpc, fabs(d_1)) < (d_2 = *stp - stpq, fabs(d_2)))
	{
	  stpf = stpc;
	} else {
	  stpf = stpq;
	}
    } else {
      if ((d_1 = *stp - stpc, fabs(d_1)) > (d_2 = *stp - stpq, fabs(d_2)))
	{
	  stpf = stpc;
	} else {
	  stpf = stpq;
	}
    }
	
    /* fourth case. a lower function value, derivatives of the */
    /* same sign, and the magnitude of the derivative does */
    /* not decrease. if the minimum is not bracketed, the step */
    /* is either stpmin or stpmax, else the cubic step is taken. */
	
  } else {
    *info = 4;
    bound = 0;
    if (*brackt) {
      theta = (fp - *fy) * 3 / (*sty - *stp) + *dy + dp;
      /* computing max */
      d_1 = fabs(theta), d_2 = fabs(*dy), d_1 = nt_max(d_2,d_1), d_2 = fabs(dp);
      s = nt_max(d_2,d_1);
      /* computing 2nd power */
      d_1 = theta / s;
      gamma = s * sqrt(d_1 * d_1 - *dy / s * (dp / s));
      if (*stp > *sty) {
	gamma = -gamma;
      }
      p = gamma - dp + theta;
      q = gamma - dp + gamma + *dy;
      r = p / q;
      stpc = *stp + r * (*sty - *stp);
      stpf = stpc;
    } else if (*stp > *stx) {
      stpf = stpmax;
    } else {
      stpf = stpmin;
    }
  }
    
  /* update the interval of uncertainty. this update does not */
  /* depend on the new step or the case analysis above. */
    
  if (fp > *fx) {
    *sty = *stp;
    *fy = fp;
    *dy = dp;
  } else {
    if (sgnd < 0.) {
      *sty = *stx;
      *fy = *fx;
      *dy = *dx;
    }
    *stx = *stp;
    *fx = fp;
    *dx = dp;
  }
    
  /*  compute the new step and safeguard it. */
    
  stpf = nt_min(stpmax,stpf);
  stpf = nt_max(stpmin,stpf);
  *stp = stpf;
  if (*brackt && bound) {
    if (*sty > *stx) {
      /* computing max */
      d_1 = *stx + (*sty - *stx) * (float).66;
      *stp = nt_min(*stp,d_1);
    } else {
      /* computing max */
      d_1 = *stx + (*sty - *stx) * (float).66;
      *stp = nt_max(*stp,d_1);
    }
  }
  return 0;
    
}

/**
 * Find a step which satisfies a sufficient decrease condition and a
 * curvature condition.
 *
 *   at each stage the subroutine updates an interval of 
 *   uncertainty with endpoints stx and sty. the interval of 
 *   uncertainty is initially chosen so that it contains a 
 *   minimizer of the modified function 
 *        f(x+stp*s) - f(x) - ftol*stp*(gradf(x)'s). 
 *
 *   if a step is obtained for which the modified function 
 *   has a nonpositive function value and nonnegative derivative, 
 *   then the interval of uncertainty is chosen so that it 
 *   contains a minimizer of f(x+stp*s). 
 *   the algorithm is designed to find a step which satisfies 
 *   the sufficient decrease condition 
 *         f(x+stp*s) .le. f(x) + ftol*stp*(gradf(x)'s), 
 *   and the curvature condition 
 *         abs(gradf(x+stp*s)'s)) .le. gtol*abs(gradf(x)'s). 
 *   if ftol is less than gtol and if, for example, the function 
 *   is bounded below, then there is always a step which satisfies 
 *   both conditions. if no step can be found which satisfies both 
 *   conditions, then the algorithm usually stops when rounding 
 *   errors prevent further progress. in this case stp only 
 *   satisfies the sufficient decrease condition. 
 *
 * @param n is a positive int input variable set to the number 
 *   of variables.
 *     fvalue is the value of f at x.
 *     x is an array of length n. on input it must contain the 
 *       base point for the line search. on output it contains 
 *       x + stp*s. 
 *     f is the function to evaluate
 *     fg (optional) is the function to evaluate with its gradient
 *     f is a variable. on input it must contain the value of f 
 *       at x. on output it contains the value of f at x + stp*s. 
 *     g is an array of length n. on input it must contain the 
 *       gradient of f at x. on output it contains the gradient 
 *       of f at x + stp*s. 
 *     s is an input array of length n which specifies the 
 *       search direction. 
 *     stp is a nonnegative variable. on input stp contains an 
 *       initial estimate of a satisfactory step. on output 
 *       stp contains the final estimate. 
 *     ftol and gtol are nonnegative input variables. 
 *       termination occurs when the sufficient decrease 
 *       condition and the directional derivative condition are 
 *       satisfied. 
 *       ftol should be smaller than 5.e-1
 *       suggested value = 1.e-4 for newton methods
 *                       = 1.e-1 for more exact line searches
 *
 *       gtol should be greater than 1.e-4 
 *       suggested value = 0.9
 *     xtol is a nonnegative input variable. termination occurs 
 *       when the relative width of the interval of uncertainty 
 *       is at most xtol. 
 *       suggested value = 2.2e-16
 *     stpmin and stpmax are nonnegative input variables which 
 *       specify lower and upper bounds for the step. (in this reverse 
 *       communication implementatin they are defined in a common 
 *       statement). 
 *       stpmin suggested value = 1.e-9
 *       stpmax suggested value = 1.e3
 *     info is an integer output variable set as follows:
 *       info = NLP_LINESEARCH_DESCENT_DIRECTION
 *         the search direction isn't a descent direction 
 *       info =-1 improper input parameters. 
 *       info = 1  the sufficient decrease condition and the 
 *                 directional derivative condition hold. 
 *       info = 2  relative width of the interval of uncertainty 
 *                 is at most xtol. 
 *       info = 4  the step is at the lower bound stpmin. 
 *       info = 5  the step is at the upper bound stpmax. 
 *       info = 6  rounding errors prevent further progress. 
 *                 there may not be a step which satisfies the 
 *                 sufficient decrease and curvature conditions. 
 *                 tolerances may be too small.
 *     work is a working vector of minimum size n.
 * 
 *   Initial version:
 *   Argonne National Laboratory. Minpack project. June 1983 
 *   Jorge J. More', David J. Thuente 
 */
int dcsrch(int n, double *grad, double *fvalue, double *x,
	   nlp_wrapper_fct wf, C_GENERIC f, int f_type,
	   double (*fg)(double *, int, void *, double *),
	   void *data, deriv1 d1,
	   double *scaling, double *s,
	   double *stp, double ftol, double xtol, 
	   double gtol, double stpmax, double stpmin,
	   double *work)
{
    /* initialized data */
  static double p66 = .66;
  static double xtrapf = 4.;
  
  /* local variables */
  double dgxm, dgym;
  int infoc;
  int info;
  double finit, width, stmin, stmax;
  int stage1;
  double width1, ftest1, dg, fm, fx, fy;
  int brackt;
  double dginit, dgtest;
  double  dgm, dgx, dgy, fxm, fym, stx, sty;

  double *newgrad;
  double suffval = 0.0, suffstep = 0.0;

  int    siter;
  int    maxiter = 10;

  newgrad = _mallocs(n, double);

  infoc = 1;
  
  /* Check the input parameters for errors. */
  if (n <= 0 || *stp <= 0.0 || ftol < 0.0 || gtol < 0.0 ||
      xtol < 0.0 || stpmin < 0.0 || stpmax < stpmin ) {
    free(newgrad);
    return NLP_LINESEARCH_IMPROPER_INPUT;
  }
  
  /* Compute the initial gradient in the search direction */
  /* and check that s is a descent direction. */
  
  dginit = cblas_ddot(n, grad, 1, s, 1);

  if (dginit >= 0.0)
    printf("\ndcsrch: Initial search direction not a descent direction\n");

  /* Initialize local variables. */
  
  brackt = 0;
  stage1 = 1;
  finit = *fvalue;
  dgtest = ftol * dginit;
  width = stpmax - stpmin;
  width1 = width / 0.5;

  /* The variables stx, fx, dgx contain the values of the step,
     function, and directional derivative at the best step.
     the variables sty, fy, dgy contain the value of the step,
     function, and derivative at the other endpoint of
     the interval of uncertainty.
     The variables stp, f, dg contain the values of the step,
     function, and derivative at the current step. */

  stx = 0.0;
  fx = finit;
  dgx = dginit;
  sty = 0.0;
  fy = finit;
  dgy = dginit;
  
  /* start of iteration. */
  
  for (siter=0; siter < maxiter; siter++) {
    /* for (siter=1; siter < maxiter; siter++) { */

    /* Set the minimum and maximum steps to correspond
       to the present interval of uncertainty. */

    if (brackt) {
      stmin = nt_min(stx,sty);
      stmax = nt_max(stx,sty);
    } else {
      stmin = stx;
      stmax = *stp + xtrapf * (*stp - stx); 
    }
    
    /* Force the step to be within the bounds stpmax and stpmin. */
    
    *stp = nt_max(*stp,stpmin);
    *stp = nt_min(*stp,stpmax);
         
    /* If an unusual termination is to occur then let
       tp be the lowest point obtained so far. */
    if ( ( brackt && (*stp <= stmin || *stp >= stmax) ) || 
	 infoc == 0 ||
	 ( brackt && stmax - stmin <= xtol * stmax) ) {
      *stp = stx;
    }

    /* Evaluate the function and gradient at stp
       and compute the directional derivative. */
    cblas_dcopy(n, x, 1, work, 1);
    cblas_daxpy(n,*stp,s,1,work,1);
    if (fg)
      *fvalue = fg(work, n, data, newgrad);
    else {
      *fvalue = wf(f, f_type, work, n, 0, data, NULL, newgrad,
		   NULL, NULL, 0, &info);
      if (d1)
	d1(work, n, data, newgrad);
      else {
	/* c_error */
	fprintf(stderr, "** Error in dcsrch routine: no derivates **\n");
	exit(1);
      }
    }

    info = 0;
    dg = cblas_ddot(n, newgrad, 1, s, 1);
    ftest1 = finit + *stp * dgtest;
 
    /*    test for convergence. */
    if ( ( brackt && (*stp <= stmin || *stp >= stmax) ) || infoc == 0)
      info = NLP_LINESEARCH_ROUNDING_ERRORS;
    else if (*stp == stpmax && *fvalue <= ftest1 && dg <= dgtest)
      info = NLP_LINESEARCH_LARGE_STEP;
    else if (*stp == stpmin && (*fvalue > ftest1 || dg >= dgtest))
      info = NLP_LINESEARCH_SMALL_STEP;
    else if (brackt && stmax - stmin <= xtol * stmax)
      info = NLP_LINESEARCH_SMALL_INTERVAL;
    else if (*fvalue <= ftest1 && fabs(dg) <= gtol * (-dginit))
      info = NLP_LINESEARCH_WOLFE_CONDITIONS;
    else if (*fvalue < ftest1) {
      cblas_dcopy(n, newgrad, 1, grad, 1);
      suffstep = *stp;
      suffval = *fvalue;
    }
    
    /*    check for termination. */
    if (info != 0) {
      /* if (suff > 0.0) *stp = suff; */
      cblas_dcopy(n, work, 1, x, 1);
      cblas_dcopy(n, newgrad, 1, grad, 1);
      free(newgrad);
      return info;
    }
    /* in the first stage we seek a step for which the modified
       function has a nonpositive value and nonnegative derivative.
       if (stage1 && *fvalue <= ftest1 && dg >=_min(ftol,gtol) * dginit) */
    if (stage1 && *fvalue <= ftest1 && dg >= 0) stage1 = 0;

    /* a modified function is used to predict the step only if
       we have not obtained a step for which the modified
       function has a nonpositive function value and nonnegative
       derivative, and if a lower function value has been
       obtained but the decrease is not sufficient. */

    if (stage1 && *fvalue <= fx && *fvalue > ftest1) {
      
      /* define the modified function and derivative values. */

      fm = *fvalue - *stp * dgtest; 
      fxm = fx - stx * dgtest;
      fym = fy - sty * dgtest;
      dgm = dg - dgtest;
      dgxm = dgx - dgtest;
      dgym = dgy - dgtest;

      /* call cstep to update the interval of uncertainty
	 and to compute the new step. */
      
      dcstep(&stx, &fxm, &dgxm, &sty, &fym, &dgym, stp, fm,
	     dgm, &brackt, stmin, stmax, &infoc); 

      /* reset the function and gradient values for f. */

      fx = fxm + stx * dgtest;
      fy = fym + sty * dgtest;
      dgx = dgxm + dgtest;
      dgy = dgym + dgtest;
    } else {
      /* call dcstep to update the interval of uncertainty
	 and to compute the new step. */

      dcstep(&stx, &fx, &dgx, &sty, &fy, &dgy, stp, *fvalue,
	     dg, &brackt, stmin, stmax, &infoc);
    }

    /*  force a sufficient decrease in the size of the
	interval of uncertainty. */

    if (brackt) {
      if (fabs(sty - stx) >= p66 * width1)
	*stp = stx + 0.5*(sty - stx);
      width1 = width;
      width = fabs(sty - stx);
    }

    /*  end of iteration. */
  }

  /* We have exhausted the maximum number of iterations */
  free(newgrad);

  /* If a sufficient decrease has been found, use it */
  if (suffstep > 0.0) {
    *stp = suffstep;
    *fvalue = suffval;
    cblas_dcopy(n, x, 1, work, 1);
    cblas_daxpy(n, *stp, s, 1, work, 1);
    cblas_dcopy(n, work, 1, x, 1);
    return NLP_LINESEARCH_SUFFICIENT_DECREASE;
  }
  else
    return NLP_LINESEARCH_FAILURE;
}
