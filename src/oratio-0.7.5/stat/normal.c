/**
 * normal.c
 *
 * STAT Library
 * Copyright (C) 1999-2007 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

/**
 * Normal distribution
 */

#include "oratio/ntools.h"

#include <math.h>
#include <float.h>

/**
 * Returns the value of a Gaussian probability distribution function
 * evaluated at x.
 *
 * @param x is the point where the probability distribution is evaluated.
 *
 * @return the value of the Gaussian pdf evaluated at x.
 */
double st_gaussian_pdf(double x)
{
  return (1/sqrt(2*M_PI))*exp(-0.5*x*x);
}

/**
 * Evaluate the density function of a normal function of mean mu and
 * standard deviation sigma, at x.
 *
 * @param x is the point where the probability distribution is evaluated.
 * @param mu is the mean of the distribution.
 * @param sigma is the standard deviation of the distribution.
 *
 * @return the value of the density fonction at x.
 */
double st_normal_pdf(double x, double mu, double sigma)
{
  return st_gaussian_pdf((x-mu)/sigma)/sigma;
}

/**
 * Returns the value of a Gaussian cumulative distribution function
 * evaluated at x, using Cody's algorithm.
 *
 * @param x is the upper limit of integration
 *
 * @return the value of a Gaussian distribution at x.
 */
double st_gaussian_cody_cdf(double x)
{
  /**
   * This routine is based on the function cumnor from DCDFLIB, which is a
   * slight modification of the function ANORM from
   *
   *  W.D. Cody (1993), "ALGORITHM 715: SPECFUN - A Portable FORTRAN
   *  Package of Special Function Routines and Test Drivers",
   *  ACM Transactions on Mathematical Software. 19, 22-32.
   *
   * ----- 
   *
   * Original comments:
   * ------------------
   *
   *  This function evaluates the normal distribution function:
   *
   *                    / x
   *             1      |    -t²/2
   *  P(x) = ---------  |   e      dt
   *         sqrt(2pi)  |
   *                    / -oo
   *
   *  The main computation evaluates near-minimax approximations derived
   *  from those in "Rational Chebyshev approximations for the error
   *  function" by W. J. Cody, Math. Comp., 1969, 631-637.
   *  This transportable program uses rational functions that theorically
   *  approximate the normal distribution to at least 18 significant
   *  decimal digits. The accuracy achieved depends on the arithmetic
   *  system, the compiler, the intrinsic functions, and proper selection
   *  of the machine-dependant constants.
   */

  static double a[5] = {
    2.2352520354606839287e00, 1.6102823106855587881e02,
    1.0676894854603709582e03, 1.8154981253343561249e04,
    6.5682337918207449113e-2
  };
  static double b[4] = {
    4.7202581904688241870e01, 9.7609855173777669322e02,
    1.0260932208618978205e04, 4.5507789335026729956e04
  };
  static double c[9] = {
    3.9894151208813466764e-1, 8.8831497943883759412e00,
    9.3506656132177855979e01, 5.9727027639480026226e02,
    2.4945375852903726711e03, 6.8481904505362823326e03,
    1.1602651437647350124e04, 9.8427148383839780218e03,
    1.0765576773720192317e-8
  };
  static double d[8] = {
    2.2266688044328115691e01, 2.3538790178262499861e02,
    1.5193775994075548050e03, 6.4855582982667607550e03,
    1.8615571640885098091e04, 3.4900952721145977266e04,
    3.8912003286093271411e04, 1.9685429676859990727e04
  };
  static double p[6] = {
    2.1589853405795699e-1, 1.274011611602473639e-1, 2.2235277870649807e-2,
    1.421619193227893466e-3, 2.9112874951168792e-5, 2.307344176494017303e-2
  };
  static double q[5] = {
    1.28426009614491121e00, 4.68238212480865118e-1, 6.59881378689285515e-2,
    3.78239633202758244e-3, 7.29751555083966205e-5
  };
  static double sqrpi = 3.9894228040143267794e-1;
  double half = 0.5e0;
  double eps, xmin;
  double del, y, xden, xnum, xsq, temp;
  int i;

  double alpha;  /* quantile value at x. */
  double calpha; /* 1 minus alpha. */

  /* Machine-dependent constants */
  eps = DBL_EPSILON*0.5e0;
  xmin = DBL_MIN;

  y = fabs(x);
  /* Evaluation dans le cas où |x| <= 0.66291 */
  if (y <= 0.66291e0) {
    xsq = 0.0;
    if (y > eps) xsq = x*x;
    xnum = a[4]*xsq;
    xden = xsq;
    for (i = 0; i < 3; i++) {
      xnum = (xnum+a[i])*xsq;
      xden = (xden+b[i])*xsq;
    }
    alpha = x*(xnum+a[3])/(xden+b[3]);
    temp = alpha;
    alpha = half+temp;
    calpha = half-temp;
  }
  /* Evaluation dans le cas où 0.66291 <= |X| <= sqrt(32) */
  else if (y <= 5.656854248e0) {
    xnum = c[8]*y;
    xden = y;
    for (i = 0; i < 7; i++) {
      xnum = (xnum+c[i])*y;
      xden = (xden+d[i])*y;
    }
    alpha = (xnum+c[7])/(xden+d[7]);
    xsq = (double)((int)((y*1.60e0)/1.60e0)); /* truncate the result */
    del = (y-xsq)*(y+xsq);
    alpha = exp(-(xsq*xsq*half))*exp(-(del*half))*alpha;
    calpha = 1.0-alpha;
    if (x > 0.0) {
      temp = alpha;
      alpha = calpha;
      calpha = temp;
    }
  }
  /** Evaluation si |x| > sqrt(32) **/
  else {
    alpha = 0.0;
    xsq = 1.0/(x*x);
    xnum = p[5]*xsq;
    xden = xsq;
    for (i = 0; i < 4; i++) {
      xnum = (xnum+p[i])*xsq;
      xden = (xden+q[i])*xsq;
    }
    alpha = xsq*(xnum+p[4])/(xden+q[4]);
    alpha = (sqrpi-alpha)/y;
    xsq = (double)((int)((x*1.60e0)/1.60e0));
    del = (x-xsq)*(x+xsq);
    alpha = exp(-(xsq*xsq*half))*exp(-(del*half))*alpha;
    calpha = 1.0-alpha;
    if (x > 0.0) {
      temp = alpha;
      alpha = calpha;
      calpha = temp;
    }
  }

  if (alpha < xmin) alpha = 0.0e0;

  return alpha;
}

/**
 * Returns the value at x of a Gaussian complementary cumulative distribution
 * (ccdf) function. Marsaglia "fast" implementation.
 *
 * @param x is the point where the ccdf has to be evaluated.
 *
 * @return the Gaussian ccdf at x.
 */
double st_gaussian_fast_ccdf(double x)
{
  /**
   * Marsaglia Complementary Normal Distribution Function
   * cPhi(x) = integral from x to infinity of exp(-.5*t^2)/sqrt(2*pi), x<15
   * 15-digit accuracy for x<15, returns 0 for x>15.
   */
  
  long double v[]={0.,.65567954241879847154L,
		   .42136922928805447322L,.30459029871010329573L,
		   .23665238291356067062L,.19280810471531576488L,
		   .16237766089686746182L,.14010418345305024160L,
		   .12313196325793229628L,.10978728257830829123L,
		   .99028596471731921395e-1L,.90175675501064682280e-1L,
		   .82766286501369177252e-1L,.76475761016248502993e-1L,
		   .71069580538852107091e-1L,.66374235823250173591e-1L};

  long double h,a,b,z,t,sum,pwr;
  int i,j;
  if(x>15.) return (0.);
  if(x<-15.) return (1.);
  j=fabs(x)+1.;
  z=j;
  h=fabs(x)-z;
  a=v[j];
  b=z*a-1.;
  pwr=1.;
  sum=a+h*b;

  for(i=2;i<60;i+=2){
    a=(a+z*b)/i;
    b=(b+z*a)/(i+1);
    pwr=pwr*h*h;
    t=sum;
    sum=sum+pwr*(a+h*b);
    if(sum==t) break; }

  sum=sum*exp(-.5*x*x-.91893853320467274178L);
  if (x<0.) sum=1.-sum;
  return ((double) sum);
}

/**
 * Returns the value of a Gaussian cumulative distribution function
 * evaluated at x.
 *
 * @param x: upper limit of integration
 *
 * @return the Gaussian cdf at x.
 */
double st_gaussian_cdf(double x)
{
  return st_gaussian_cody_cdf(x);
}

/**
 * Returns the value at x of a Gaussian complementary cumulative distribution
 * (ccdf) function. Marsaglia "fast" implementation.
 *
 * @param x is the point where the ccdf has to be evaluated.
 *
 * @return the Gaussian ccdf at x.
 */
double st_gaussian_ccdf(double x)
{
  return st_gaussian_cody_cdf(-x);
}

/**
 * Returns the value at x of a Gaussian cumulative distribution
 * (cdf) function. Marsaglia "fast" implementation.
 *
 * @param x is the point where the cdf has to be evaluated.
 *
 * @return the Gaussian cdf at x.
 */
double st_gaussian_fast_cdf(double x)
{
  long double s=x, t=0, b=x, q=x*x, i=1;
 
  if (x < -8.0) return st_gaussian_ccdf(-x);
  else if (x > 8.3) return 1.0;  /* truncation error does not allow a
				  better precision. */

  /* We are in the "good" range. */
  while (s!=t) s=(t=s)+(b*=q/(i+=2));
  
  return .5+s*exp(-.5*q-.91893853320467274178L);
}

/**
 * Compute an approximation of the cumulative function of a normal
 * distribution, of mean mu and standard deviation sigma, at x.
 */
double st_normal_cdf(double x, double mu, double sigma)
{
  return st_gaussian_cdf((x-mu)/sigma);
}

/**
 * Compute the Gaussian inverse cumulative distribution function.
 *
 * @param u is the value corresponding to the desired quantile.
 *
 * @return the value of the inverse cumulative function of a
 *   gaussian distribution evaluated at p.
 */
double st_gaussian_icdf(const double u)
{
  /**
   * Ref.:
   *  M.J. Wichura (1988), Algorithm AS 241: The Percentage Points
   *  of the Normal Distribution, Applied Statistics, 37, 477-484.
   */

  double dp, r, x;

  static const double a[8] = {
    3.3871328727963666080E0,
    1.3314166789178437745E2,
    1.9715909503065514427E3,
    1.3731693765509461125E4,
    4.5921953931549871457E4,
    6.7265770927008700853E4,
    3.3430575583588128105E4,
    2.5090809287301226727E3
  } ;

  static const double b[7] = {
    4.2313330701600911252E1,
    6.8718700749205790830E2,
    5.3941960214247511077E3,
    2.1213794301586595867E4,
    3.9307895800092710610E4,
    2.8729085735721942674E4,
    5.2264952788528545610E3
  };

  static const double c[8] = {
    1.42343711074968357734E0,
    4.63033784615654529590E0,
    5.76949722146069140550E0,
    3.64784832476320460504E0,
    1.27045825245236838258E0,
    2.41780725177450611770E-1,
    2.27238449892691845833E-2,
    7.74545014278341407640E-4
  };

  static const double d[7] = {
    2.05319162663775882187E0,
    1.67638483018380384940E0,
    6.89767334985100004550E-1,
    1.48103976427480074590E-1,
    1.51986665636164571966E-2,
    5.47593808499534494600E-4,
    1.05075007164441684324E-9
  };

  static double e[8] = {
    6.65790464350110377720E0,
    5.46378491116411436990E0,
    1.78482653991729133580E0,
    2.96560571828504891230E-1,
    2.65321895265761230930E-2,
    1.24266094738807843860E-3,
    2.71155556874348757815E-5,
    2.01033439929228813265E-7
  };

  static const double f[7] = {
    5.99832206555887937690E-1,
    1.36929880922735805310E-1,
    1.48753612908506148525E-2,
    7.86869131145613259100E-4,
    1.84631831751005468180E-5,
    1.42151175831644588870E-7,
    2.04426310338993978564E-15
  };

  /* Produce infinite values if u is set to 0.0 or 1.0 */
  if (u == 1.0) return nt_infinity();
  else if (u == 0.0) return -nt_infinity();

  dp = u-0.5;

  if (fabs(dp) < 0.425) {
    /* small values */
    r = 0.180625-dp*dp;
    return dp*(((((((a[7]*r + a[6])*r + a[5])*r + a[4])*r + a[3])*r
		+ a[2])*r + a[1])*r + a[0]) /
      (((((((b[6]*r + b[5])*r + b[4])*r + b[3])*r + b[2])*r + b[1])*r
	+ b[0])*r + 1.0);     
  }

  r = (dp < 0.0) ? u : 1.0 - u;
  r = sqrt(-log(r));

  if (r <= 5.0) {
    /* intermediate values */
    r -= 1.6;
    x = (((((((c[7]*r + c[6])*r + c[5])*r + c[4])*r
	    + c[3])*r + c[2])*r + c[1])*r + c[0]) /
      (((((((d[6]*r + d[5])*r + d[4])*r + d[3])*r
	  + d[2])*r + d[1])*r + d[0])*r + 1.0);
  }
  else {
    /* tail values */
    r -= 5.0;
    x = (((((((e[7]*r + e[6])*r + e[5])*r + e[4])*r + e[3])*r + e[2])*r
	  + e[1])*r + e[0]) /
      (((((((f[6]*r + f[5])*r + f[4])*r + f[3])*r
	  + f[2])*r + f[1])*r + f[0])*r + 1.0);
  }

  if (dp < 0.0)
    return -x;
  else
    return x;
}

/**
 * Compute the Normal inverse cumulative distribution function.
 *
 * @param p is the value (0 <= p <= 1) with respect we want to compute
 *   the associated quantile.
 * @param mu is the mean of the normal distribution.
 * @param sigma is the standard deviation of the normal distribution.
 *
 * @return the value of the inverse cumulative function of a
 *   normal distribution of mean mu and standard deviation sigma.
 */
double st_normal_icdf(double p, double mu, double sigma)
{
  return mu + sigma*st_gaussian_icdf(p);
}

/**
 * Return the mean of a normal distribution.
 *
 * @param mu is the mean of the normal distribution.
 * @param sigma is the standard deviation of the normal distribution.
 *
 * @return the mean of the normal distribution.
 */
double st_normal_mean(double mu, double sigma)
{
  return mu;
}
