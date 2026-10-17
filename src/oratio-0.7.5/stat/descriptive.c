/**
 * descriptive.c
 *
 * ORATIO Library
 * Copyright (C) 1999-2013 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include <math.h>

#include <oratio/ntools.h>

/**
 * Update the mean with one new observation.
 *
 * @param n the number of observations used to compute the previous mean.
 * @param mean is the previous mean, computed without the current observation.
 * @param obs is the current observation.
 *
 * @return the value of the updated mean, combining the previous mean and
 *   the current observation.
 */
double st_mean_update(int n, double mean, double obs)
{
  return mean+(obs-mean)/(n+1);
}

/**
 * Remove one observation from the mean.
 *
 * @param n the number of observations used to compute the previous mean.
 * @param mean is the previous mean, computed without the current observation.
 * @param obs is the current observation.
 *
 * @return the value of the updated mean, combining the previous mean and
 *   the current observation.
 */
double st_mean_remove(int n, double mean, double obs)
{
  if (n <= 1) return nt_nan();

  return mean+(mean-obs)/(n-1);
}

/**
 * Compute the mean of observations obs[l*inc] to obs[(u-1)*inc].
 * If u <= l, return 0.
 *
 * @param obs is a pointer to an array of sorted observations.
 * @param inc is the increment to use in the obs array (usually 1).
 * @param l is the index of the lowest observation to take into
 *   consideration.
 * @param u is the index of the first observation to discard from
 *   consideration.
 *
 * @return 0 if u <= l, return 0, otherwise the mean of
 *   observations obs[l*inc] to obs[(u-1)*inc] (included).
 */
double st_mean(double *obs, int inc, int l, int u)
{
  int i;
  double mean = 0.0;

  for (i = l; i < u; i++) mean += (obs[i*inc]-mean)/(i-l+1);

  return mean;
}

/**
 * Compute the median from a sorted observations array.
 *
 * @param obs is a pointer to an array of sorted observations.
 * @param inc is the increment to use in the obs array (usually 1).
 * @param low is the index of the lowest observation to take into
 *   consideration.
 * @param up is the index of the first observation to discard from
 *   consideration. 
 *
 * @return the median of [ obs[low], obs[up] [. If the number of
 *   observations is even, return the average between the last element
 *   in the lower half part and the first element in the upper part of
 *   the observations, as usually prescribed.
 */
double st_median(double *obs, int inc, int low, int up)
{
  int n = up-low;
  int pos = n >> 1;
  double x;

  if (n % 2 == 0)
    /* n is even */
    x = (obs[pos]+obs[(pos - 1)*inc])*0.5;
  else
    /* n is odd */
    x = obs[pos];

  return x;
}

/**
 * Compute the median from a sorted observations array, restricting
 * the median to be one of the present observations.
 *
 * @param obs is a pointer to an array of sorted observations.
 * @param inc is the increment to use in the obs array (usually 1).
 * @param l is the index of the lowest observation to take into
 *   consideration.
 * @param u is the index of the first observation to discard from
 *   consideration. 
 *
 * @return the median of [ obs[low], obs[up] [. If the number of
 *   elements is even, the function returns the first element in the
 *   upper half part of the observations.
 */
double st_median_observation(double *obs, int inc, int l, int u)
{  
  return obs[(int)((rint((u-l)*0.5))*inc)];
}

/**
 * Compute the variance of observations obs[low] to obs[up-1].
 *
 * @param mean is the mean of observations obs[low] to obs[up-1]
 *   (included).
 * @param obs is a pointer to an array of observations.
 * @param inc is the increment to use in the obs array (usually 1).
 * @param low is the index of the lowest observation to take into
 *   consideration.
 * @param up is the index of the first observation to discard from
 *   consideration. 
 *
 * @return 0 is idx_high <= idx_low+1 (i.e. there strictly less than one
 *   observation), otherwise the variance of observations obs[low] to
 *   obs[up-1] (included).
 */
double st_variance(double mean, double *obs, int inc, int low, int up)
{
  int i;
  int n = up-low;;
  double variance = 0.0, delta;

  if (up > low+1) {
    for (i = low; i < up; i++) {
      delta = obs[i*inc]-mean;
      variance += (delta*delta-variance)/(i-low+1);
    }
  }

  return variance*((double)n/(double)(n-1));
}

/**
 * Compute the mean and variance of the observations going from
 * obs[low] to obs[up-1].
 *
 * @param obs is a pointer to an array of observations.
 * @param inc is the increment to use in the obs array (usually 1).
 * @param l is the index of the lowest observation to take into
 *   consideration.
 * @param u is the index of the first observation to discard from
 *   consideration. 
 */
void st_mean_variance(double *obs, int inc, int l, int u,
			     double *mean, double *variance)
{
  *mean = st_mean(obs, inc, l, u);
  *variance = st_variance(*mean, obs, inc, l, u);
}

/**
 * Compute a quantile from empirical data.
 *
 * @param n is the dimension of the array values.
 * @param x is the (sorted) array of observations.
 * @param incx is the increment to use in the array x (usually 1).
 * @param quantile is the quantile to evaluate.
 *   to the quantiles values on output.
 *
 * @return the quantile value.
 */
double st_quantile(int n, double *x, int incx,
			  double quantile)
{
  return x[((int)(ceil(n*quantile)-1))*incx];
}

/**
 * Compute the quantiles from empirical data.
 *
 * @param n is the dimension of the array values.
 * @param x is the (sorted) array of observations.
 * @param inc is the increment to use in the array x (usually 1).
 * @param m is the number of quantiles to computes. Quantiles are
 *   considered in ]0,1[.
 * @param quantiles is an double array of m elements, whose content is
 *   set on output to the quantile values.
 */
void st_quantiles(int n, double *x, int inc,
		  int m, double *quantiles)
{
  int i;
  int idx;

  for (i = 0; i < m; i++) {
    idx = ceil(n*(i+1.0)/(m+1.0))-1;
    quantiles[i] = x[idx*inc];
  }
}

/**
 * Compute the covariance between observations vectors x and y (expert
 * routine).
 *
 * @param n is the number of observations in x and y.
 * @param x is the first observations vector.
 * @param incx is the increment to use in the array x.
 * @param meanx is the mean of vector x.
 * @param y is the second observations vector.
 * @param incy is the increment to use in the array y.
 * @param meany is the mean of vector y.
 * @param l is the index of the lowest observation to take into
 *   consideration.
 * @param u is the index of the first observation to discard from
 *   consideration. 
 *
 * @return the covariance between x and y.
 */
double st_covariance(int n, double *x, int incx, double meanx,
		     double *y, int incy, double meany,
		     int l, int u)
{
  int i;
  double covariance = 0.0;

  if (u > l+1)
    for (i = l; i < u; i++)
      covariance += (((x[i*incx]-meanx)*(y[i*incy]-meany))-covariance)/(i-l+1);

  return ((double)n/(double)(n-1))*covariance;
}
