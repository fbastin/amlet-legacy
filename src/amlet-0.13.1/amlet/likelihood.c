/**
 * likelihood.c
 *
 * AMLET Library
 * Copyright (C) 2002-2011 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include <math.h>
#include <float.h>
#include <stdio.h>

#include <amlet/amlet_config.h>
#include <amlet/dcl.h>

#include <ophelia/bsplines.h>
#include <ophelia/nonparametric.h>

#include <oratio/oratio.h>

#include "logit.h"
#include "mixedlogit.h"
#include "utilities.h"

#define DFLT_ERR   1e100  /* Default error (a value that is assumed to be much
			     larger than estimated accuracy for practical
			     applications */
#define MAX_VALUE  1e25   /* Too small values lead to bad behaviour of the
			     BFGS, a well as too large :( */
#define SAFE_VALUE 1e-150 /* this value will be used as the minimum
			     value in some computations. */

void set_L(int n, double *L, double *x)
{
  int i, j, k;

  k = 0;
  for (i = 0; i < n; i++)
    for (j = 0; j <=i; j++)
      L[i*n+j] = x[k++];
}

/**
 * @param i is the individual index.
 */
void derive_logit(int i, int n, double ***obs, int choice,
			 double den, int *avail, double *tgrad, double *grad,
			 dcl_data *md, double *utility, int shift)
{
  double prob = 0.0;
  int j;

  for (j = 0; j < n; j++) tgrad[j] = 0.0;
  /* Compute the gradient contribution from user's choice. */
  derive_utility(n, i, choice, tgrad+n, obs);
  /* Compute the gradient contribution of the other alternatives. */
  for (j = 0; j < md->nalt; j++) {
    prob = utility[j]/den;
    
    if (j != choice && avail[j]) {
      derive_utility(n, i, j, tgrad+2*n, obs);
      cblas_daxpy(n, prob, tgrad+2*n, 1, tgrad, 1);
      cblas_daxpy(n, -prob, tgrad+n, 1, tgrad, 1);
    }
  }

  for (j = 0; j < n; j++) {
    md->sgrad[i][j] = tgrad[j];
    grad[j] += (tgrad[j]-grad[j]) / (i+1-shift);
  }
}

/**
 * Set the parameters vector beta.
 *
 * @param indidx is the individual index.
 * @param s is the sampling index.
 * @param beta is on output the coefficients vector of the utility function.
 * @param theta is the vector of model parameters.
 * @param antithetic is equal to 1 if we work with antithetic variates, 0 else.
 */
void set_beta(int nvar, double *beta, double *theta,
		     const int *btype,
		     double **sample, int s, int indidx, int ssamp,
		     int size_elem, dcl_data *md,
		     int antithetic)
{
  int i, j = 0;
  int pos = 0, spos = 0;
  double draw;

  DS_List *np_list = md->np_list;
  void *obj;

  for (i = 0; i < nvar; i++) {
    if (amlet_is_nonparametric(btype[i])) {
      obj = np_list->data;
      np_list = np_list->next;
    }
    else
      obj = NULL;
    if (btype[i] == FIXED)
      beta[i] = md->fixed[j++];
    else {
      draw = 0.0;
      if (btype[i] != CONSTANT)
	draw = sample[indidx][s*md->n_random_coeff+spos++];
	if (antithetic) draw = 1-draw;
    pos += moca_set_realisation(beta+i, theta, draw, btype[i],
				  pos, obj);
    }
  }
}

void set_beta_mvn(int n, double *beta, double *theta,
			 const int *btype,
			 double **sample, int s, int indidx, int ssamp,
			 int size_elem, dcl_data *md,
			 int antithetic)
{
  double L[n*n];

  set_L(n, L, theta);
  /*  for (i = 0; i < n; i++)
      beta[i] = st_gaussian_icdf(sample[indidx][s*n+i]); */

  st_mvn_convert(CblasRowMajor, n, L, n, theta+(n*(n+1))/2, 1, beta, 1);
}

void set_beta_copula(int n, double *beta, double *theta,
			    const int *btype,
			    double **sample, int s, int indidx, int ssamp,
			    int size_elem, dcl_data *md,
			    int antithetic)
{
  double x[n], L[n*n], u[n];
  int i;

  DS_List *np_list = md->np_list;
  void *obj;

  int pos = (n*(n+1))/2;

  set_L(n, L, theta);

  for (i = 0; i < n; i++) u[i] = 0.0;
  for (i = 0; i < n; i++) x[i] = sample[indidx][s*n+i];
  //  for (i = 0; i < n; i++) x[i] = st_gaussian_icdf(sample[indidx][s*n+i]);
  st_mvn_convert(CblasRowMajor, n, L, n, u, 1, x, 1);

  for (i = 0; i < n; i++) u[i] = st_gaussian_cdf(x[i]);

  for (i = 0; i < n; i++) {
    if (amlet_is_nonparametric(btype[i])) {
      obj = np_list->data;
      np_list = np_list->next;
    }
    else
      obj = NULL;
    pos += moca_set_realisation(beta+i, theta, st_gaussian_icdf(u[i]),
				btype[i], pos, obj);
  }
}

/**
 * DERIVATIVES
 */

/**
 * Compute the derivative associated to a coefficient inside a linear
 * utility.
 *
 * @param factor is the factor in the (linear) utility expression, related to
 *   the parameter to derive.
 * @param np is a pointer to some nonparametric structure,
 *   corresponding to type, or NULL is no such object is needed.
 */
void deriv_coeff(int type, double *grad, double coeff,
			double factor, double draw, void *np)
{
  /* We should be able to avoid the variable pos here... */
  switch (type) {	
  case ERRC:
    *grad += factor*draw;
    break;
  case NORMAL:
  case UNIFORM:
  case TRIANGULAR:
    *grad += factor;
    grad[1] += factor*draw;
    break;
  case LOGNORMAL:
    *grad += coeff*factor;
    grad[1] += coeff*factor*draw;
    break;
  case RAN_SP_SPLINE: {
    NT_BSpline *bs = (NT_BSpline *)np;
    /*    we should word with a working vector in order to correctly
	  made the addition instead of a setting */
    amlet_bsplines_curve_point_grad(bs->m, bs->n, bs->p, draw,
				    bs->U, bs->N, grad, factor);
    /* multiply by factor */
    //    cblas_dscal(bs->m, factor, grad, 1);
    break;
  }
  case CONSTANT:
  default:
    *grad += factor;
  }
}

/**
 * Compute the utilities derivatives.
 *
 * @param ind is the individual index.
 * @param samp_idx is the sampling draw index.
 * @param choice is the index of the alternative chosen by ind.
 * @param md is the data structure about the problem.
 */
void sll_derivate_utilities(int ind, int time, int samp_idx,
				   int choice, const double *beta,
				   const double *prob,
				   double *tgrad, dcl_data *md)
{
  int l, k;
  int pos, spos;
  int obsidx;       /* observation index */

  int *avail = md->availability[md->rep[ind]+time];

  DS_List *tmp = md->np_list;
  void *obj;

  obsidx = md->rep[ind] + time;

  pos = spos = 0;

  for (l = 0; l < md->nvar; l++) {
    if (md->btype[l] != FIXED) {
      obj = NULL;

      if (amlet_is_nonparametric(md->btype[l])) {
	obj = tmp->data;
	tmp = tmp->next;
      }

      for (k = 0; k < md->nalt; k++) {
	if (k != choice && avail[k]) {
	  deriv_coeff(md->btype[l], tgrad+pos, beta[l],
		      prob[k]*
		      (md->obs[obsidx][k][l]-md->obs[obsidx][choice][l]),
		      (md->btype[l] == CONSTANT) ? 0.0 :
		      md->sampling[ind][samp_idx*md->n_random_coeff+spos],
		      obj);
	}
      }
      pos += moca_count_parameters(md->btype[l], obj);
      if (md->btype[l] != CONSTANT) spos++;
    }
  }
}

/**
 * Compute the utilities derivatives.
 *
 * @param ind is the individual index.
 * @param samp_idx is the sampling draw index.
 * @param choice is the index of the alternative chosen by ind.
 * @param md is the data structure about the problem.
 */
void sll_derivate_utilities_mvn(int ind, int time, int samp_idx,
				       int choice, const double *beta,
				       const double *prob,
				       double *tgrad, dcl_data *md)
{
  int i, l, k;
  int pos, spos, pos2, n;
  int obsidx;       /* observation index */

  int *avail = md->availability[md->rep[ind]+time];

  DS_List *tmp = md->np_list;
  void *obj;

  double draw, factor;

  obsidx = md->rep[ind] + time;

  pos = 0;
  pos2 = (n*(n+1))/2;

  for (l = 0; l < md->nvar; l++) {

    for (k = 0; k < md->nalt; k++) {
      if (k != choice && avail[k]) {

	factor = prob[k]*(md->obs[obsidx][k][l]-md->obs[obsidx][choice][l]);
	
	for (i = 0; i <= l; i++) {
	  draw = md->sampling[ind][samp_idx*md->n_random_coeff+i];
	  tgrad[pos] += draw*factor;
	}

	tgrad[pos2] += factor;
      }

      pos++;
      pos2++;
    }
  }
}

/**
 * Compute the utilities of the various alternatives.
 *
 * @param obs it a tri-dimensional array with the observations. The
 *   first index is the observation index, the secont the alternative
 *   and the last one the parameter index.
 * @param beta is the coefficients vector of the utility functions.
 * @param time is the time repetition index (0 if not applicable)
 * @param md is the data structure about the problem.
 * @param den is on output the denominator of the logit function.
 */
void compute_utilities(double *beta, double ***obs,
			      dcl_data *md, int ind, int time,
			      double *utility, double *den)
{
  int i, k;
  int *indavail;      /* alternatives availability */

  *den = 0.0;

  i = (md->n_random_coeff) ? md->rep[ind]+time : ind;
  indavail = md->availability[i];

  /* compute the utility of each alternative */
  for (k = 0; k < md->nalt; k++) {
    if (indavail[k])
      utility[k] = cblas_ddot(md->nvar, beta, 1,
			      obs[i][k], 1);
  }

  amlet_utilities_safeguard(utility, den, indavail, md->nalt, -1);
}

/**
 * Compute the probability of the observed choice (for the particular draw).
 *
 * @param theta is the vector of model parameters.
 * @param md is the data structure about the problem.
 * @param s: sampling draw index
 * @param ind: individu index
 * @param derivate: 1 if the derivatives have to be computed, 0 else
 *
 * @return the choice probability.
 */ 
double compute_prob_mixed(double *theta, double ***obs,
				 double **sampling, dcl_data *md,
				 int s, int ind, int max_size,
				 double *utility, int derivate, int antithetic)
{
  double den;
  double *single_prob;           /* probabilities for one repetition */
  double prob = 1.0;             /* probability of choice over the
				    repetitions */
  double *beta;
  double *tgrad;

  int time = 0;                  /* repetition index */
  int k;                         /* loop index */

  int size;

  size = (sizeof(double) == sizeof(float)) ? FLOATBITS : DOUBLEBITS;

  single_prob = _mallocs(md->nalt, double);
  beta = _mallocs(md->nvar, double);
  tgrad = _mallocs(md->npar, double);

  set_beta(md->nvar, beta, theta, md->btype, sampling,
	   s, ind, max_size, size, md, antithetic);

  /* initialize the partial gradient vector to zero */
  for (k = 0; k < md->npar; k++) tgrad[k] = 0.0;

  for (time = 0; time < md->rep[ind+1]-md->rep[ind]; time++) {
    compute_utilities(beta, obs, md, ind, time, utility, &den);

    for (k = 0; k < md->nalt; k++)
      single_prob[k] = utility[k]/den;

    if (derivate)
      sll_derivate_utilities(ind, time, s, md->choices[md->rep[ind]+time],
			     beta, single_prob, tgrad, md);

    prob *= single_prob[md->choices[md->rep[ind]+time]];
  }

  /* update md->sgrad[ind] */
  if (derivate) {
    for (k = 0; k < md->npar; k++)
      md->sgrad[ind][k] += (prob*tgrad[k]-md->sgrad[ind][k])/(s+1);
    if (md->bias) {
      for (k = 0; k < md->npar; k++)
	md->somgrad[ind][k] +=
	  (2*prob*prob*tgrad[k] - md->somgrad[ind][k])/(s+1);
    }
  }

  free(single_prob);
  free(beta);
  free(tgrad);

  return prob;
}

/**
 * Compute the opposite of the log-likelihood and its gradient.
 *
 * @param theta is the vector of model parameters.
 * @param md is the data structure about the problem.
 *
 * @return 0.0.
 */
double likelihood_gradient(double *theta, int n, void *data, double *grad)
{
  int i, j;
  double den;
  dcl_data *md = (dcl_data *)data;

  double *utility, *tgrad;

  utility = _mallocs(md->nalt, double);
  tgrad = _mallocs(3*n, double);

  for (i = 0; i < n; i++) grad[i] = 0.0;

  for (i = 0; i < md->nobs; i++) {
    for (j = 0; j < md->nalt; j++) utility[j] = 0.0;

    /* compute the utility of each alternative */
    compute_utilities(theta, md->obs, md, i, 0, utility, &den);

    derive_logit(i, n, md->obs, md->choices[i], den, md->availability[i],
		 tgrad, grad, md, utility, 0);
  }

  free(utility);
  free(tgrad);

  return 0.0;
}

/**
 * Multinomial logit log-likelihood functions
 */

/**
 * Compute the opposite of the log-likelihood function
 *
 * @param theta is the vector of model parameters.
 * @param choices is a constant integer array describing the
 *   individual choices. More precisely, choices[i] is the index of
 *   the alternative chosen by individual i.
 * @param avail is a doubly indexed integer array indicating the
 *   individual alternatives availability: avail[i][j] is equal to 1 if
 *   the jth alternative is available for individual i, 0 else.
 * @param md is the data structure about the problem.
 * @param shift designs the shift in the individuals count.
 */
double likelihood(double *theta, int np,
		  double ***obs, const int *choices,
		  int **avail, int nvar, dcl_data *md,
		  int shift)
{

  int i, j;
  double ll, prob, den;

  double *utility;

  utility = _mallocs(md->nalt, double);

  ll = 0.0;

  for (i = shift; i < np+shift; i++) {
    for (j = 0; j < md->nalt; j++) utility[j] = 0.0;

    /* compute the utility of each alternative */
    compute_utilities(theta, obs, md, i, 0, utility, &den);

    prob = utility[choices[i]]/den;

    md->sprob[i] = prob;

    ll += (log(prob)-ll)/(i+1-shift);
  }

  free(utility);

  return -ll;
}

/**
 * Compute the opposite of the log-likelihood function.
 * 
 * @param theta is the vector of model parameters.
 * @param choices is a constant integer array describing the
 *   individual choices. More precisely, choices[i] is the index of
 *   the alternative chosen by individual i.
 * @param avail is a doubly indexed integer array indicating the
 *   individual alternatives availability: avail[i][j] is equal to 1 if
 *   the jth alternative is available for individual i, 0 else.
 * @param md is the data structure about the problem.
 */
double lg(double *theta, int np,
	  double ***obs, const int *choices,
	  int **avail, int nvar, double *grad, dcl_data *md,
	  int shift)
{
  int i, j;
  double ll, prob, den;
  double *utility, *tgrad;

  utility = _mallocs(md->nalt, double);
  tgrad = _mallocs(3*nvar, double);

  ll = 0.0;

  for (i = 0; i < nvar; i++) grad[i] = 0.0;

  for (i = shift; i < np+shift; i++) {
    for (j = 0; j < md->nalt; j++) utility[j] = 0.0;

    /* compute the utility of each alternative */
    compute_utilities(theta, obs, md, i, 0, utility, &den);

    prob = nt_max(utility[choices[i]]/den, SAFE_VALUE);
    md->sprob[i] = prob;
    ll += (log(prob)-ll)/(i+1-shift);

    derive_logit(i, nvar, obs, choices[i], den, avail[i],
		 tgrad, grad, md, utility, shift);
  }

  free(utility);
  free(tgrad);

  return -ll;
}

/**
 *
 * @param md is the data structure about the problem.
 */
double sll_variance_draw(int ind, int size, dcl_data *md)
{
  /* Compute sigma^2/mu^2, with sigma^2 = E[L^2]-mu*mu */
  return md->som[ind]/(md->sprob[ind]*md->sprob[ind])-1.0;
}

/**
 * SIMULATED LIKELIHOOD
 */

/**
 * Compute the opposite of the simulated log-likelihood function
 *
 * @param theta is the vector of model parameters.
 * @param choices is a constant integer array describing the
 *   individual choices. More precisely, choices[i] is the index of
 *   the alternative chosen by individual i.
 * @param md is the data structure about the problem.
 */
double simulated_likelihood(double *theta, double ***obs,
			    const int *choices, double **sampling, int nvar,
			    int ssamp, int max_size, dcl_data *md)
{
  int i, j, k;
  double ll, prob, prob_ind;
  double *utility;

  double variance = 0.0, err = 0.0;

  utility = _mallocs(md->nalt, double);

  ll = 0.0;

  for (i = 0; i < md->np; i++) {
    prob = 0.0;
    md->sprob[i] = 0.0;
    md->som[i] = 0.0;  

    for (k = 0; k < md->nalt; k++) utility[k] = 0.0;

    for (j = 0; j < ssamp; j++) {
      prob_ind = compute_prob_mixed(theta, obs, sampling, md, j, i,
				    max_size, utility, 0, 0);
      /* add to the simulated probability */
      md->sprob[i] += (prob_ind - md->sprob[i]) / (j+1);
      /* update the second order moment */
      md->som[i] += (prob_ind*prob_ind-md->som[i]) / (j+1);
      prob += (prob_ind-prob)/(j+1);
    }

    /* we check the value of prob to avoid failure of the
       logarithm operation */
    prob = nt_max(prob, SAFE_VALUE);

    ll += (log(prob)-ll)/(i+1);
    /* compute the variance */
    variance += sll_variance_draw(i, ssamp, md);
  }

  variance /= ssamp-1;
  err = md->err_quantile/md->np*sqrt(variance);
  /* bias correction */
  if (md->bias) ll -=  simulation_bias(md->np, md->err_quantile, err);

  free(utility);

  md->neval += ssamp;

  return -ll;
}

/**
 * Compute the opposite of the log-likelihood function.
 *
 * @param theta is the vector of model parameters.
 * @param choices is a constant integer array describing the
 *   individual choices. More precisely, choices[i] is the index of
 *   the alternative chosen by individual i.
 * @param md is the data structure about the problem.
 * @param err is a pointer to a double whose content is undefined on
 *   input, and corresponds to the estimated accuracy on output, at
 *   the level md->level.
 *
 * Compute the opposite of the log-likelihood function.
 */
double slg(double *theta, double ***obs, const int *choices,
	   double **sampling, const int *btype, int ssamp,
	   int max_size, dcl_data *md, double *err, double *grad)
{
  int i, j, k, m;
  double ll;

  double *utility, prob;

  double variance = 0.0;
  double *var = NULL;

  if (md->bias) var = _mallocs(md->np, double);

  utility = _mallocs(md->nalt, double);

  for (i = 0; i < md->npar; i++) grad[i] = 0.0;

  ll = 0.0;
  m = 0;

  for (i = 0; i < md->np; i++) {

    for (k = 0; k < md->nalt; k++) utility[k] = 0.0;
    md->sprob[i] = 0.0;
    md->som[i] = 0.0;  
    for (j = 0; j < md->npar; j++) md->sgrad[i][j] = 0.0;
    if (md->bias) {
      for (j = 0; j < md->npar; j++) md->somgrad[i][j] = 0.0;
      var[i] = 0.0;
    }

    for (j = 0; j < ssamp; j++) {
      prob = compute_prob_mixed(theta, obs, sampling, md, j, i, max_size,
				utility, 1, 0);

      /* add to the simulated probability */
      md->sprob[i] += (prob - md->sprob[i]) / (j+1);
      /* update the second order moment */
      md->som[i] += (prob*prob-md->som[i]) / (j+1);
    }

    /* we check the value of prob to avoid failure of the
       logarithm operation */
    md->sprob[i] = nt_max(md->sprob[i], SAFE_VALUE);

    if (md->bias) {
      /* compute the variance */
      var[i] = sll_variance_draw(i, ssamp, md);
      variance += var[i];
    }

    /* gradient and likelihood contributions */
    m += md->rep[i+1]-md->rep[i];
    for (k = 0; k < md->npar; k++) {
      if (md->bias)
	grad[k] += (md->sgrad[i][k]/md->sprob[i] - grad[k]
		    + (md->somgrad[i][k]/(2*md->sprob[i]*md->sprob[i])
		       - md->sgrad[i][k]*(var[i]+1)/md->sprob[i])/ssamp)/(i+1);
      else
	grad[k] += (md->sgrad[i][k]/md->sprob[i] - grad[k])/(i+1);
    }

    ll += (log(md->sprob[i])-ll) / (i+1);
  }

  if (md->bias) {
    variance /= ssamp-1;
    if (err) *err = md->err_quantile/md->np*sqrt(variance);

    ll -= simulation_bias(md->np, md->err_quantile,
			  md->err_quantile/md->np*sqrt(variance));
    free(var);
  }
  else
    if (err) *err = DFLT_ERR;

  free(utility);

  return -ll;
}

/**
 * SIMULATED LIKELIHOOD ACCURACY
 */

/**
 * Evaluate the accuracy of the simulated likelihood function.
 *
 * @param theta is the vector of model parameters.
 * @param md is the data structure about the problem.
 * @param accuracy is a pointer to a double variable whose content is
 *   unknown on input, and corresponds to the accuracy of the
 *   simulated likelihood function on output.
 * @param bias is a pointer to a double variable whose content is
 *   unknown on input, and corresponds to the bias in the simulated
 *   likelihood function on output.
 */
void sll_accuracy(double *theta, dcl_data *md,
		  double *accuracy, double *bias, int start)
{
  double mean, var;

  double *slls;
  int shift;

  int i;

  switch (md->method) {
  case MOCA_MONTE_CARLO:
    /* Simulation accuracy */
    var = 0.0;
    for (i = 0; i < md->np; i++)
      var += sll_variance_draw(i, md->size_sampling, md);
    var /= md->size_sampling-1;
    *accuracy = md->err_quantile/md->np*sqrt(var);
    break;
  default:
    slls = _mallocs(md->qmcrep, double);
    shift = md->max_size*md->np;

    for (i = 0; i < md->qmcrep; i++) {
      slls[i] = simulated_likelihood(theta, md->obs, md->choices,
				     md->sampling+i*shift, md->nvar,
				     md->size_sampling, md->max_size, md);
      /* Compute the variance over the various estimations */
      st_mean_variance(slls, 1, 0, md->qmcrep, &mean, &var);
    }
    free(slls);
    /* Simulation accuracy */
    *accuracy = md->err_quantile*sqrt(var);
  }

  /* Simulation bias */
  *bias = simulation_bias(md->np, md->err_quantile, *accuracy);
}

/**
 * SIMULATED LIKELIHOOD WITH ACCURACY ESTIMATION
 */

/**
 * Compute the opposite of the log-likelihood function
 *
 * @param theta is the vector of model parameters.
 * @param choices is a constant integer array describing the
 *   individual choices. More precisely, choices[i] is the index of
 *   the alternative chosen by individual i.
 * @param md is the data structure about the problem.
 * @param err is a pointer to a double whose content is undefined on
 *   input, and corresponds to the estimated accuracy on output, at
 *   the level md->level.
 */
double sldag(double *theta, double ***obs, const int *choices,
	     double **sampling, int ssamp, int max_size,
	     dcl_data *md, double *err, double *variance,
	     double *grad, int start)
{
  int i, j, k;
  double ll;
  double *utility, *var;
  double prob;

  int *indavail;

//  FILE *f = fopen("probs.txt", "w");

  utility = _mallocs(md->nalt, double);
  var = _mallocs(md->np, double);

  ll = 0.0;
  *variance = 0.0;

  for (i = 0; i < md->npar; i++) grad[i] = 0.0;

  for (i = 0; i < md->np; i++) {
    indavail = md->availability[i];
    for (k = 0; k < md->nalt; k++) if (!indavail[k]) utility[k] = 0.0;

    if (start == 0) {
      for (j = 0; j < md->npar; j++) md->sgrad[i][j] = 0.0;
      if (md->bias) for (j = 0; j < md->npar; j++) md->somgrad[i][j] = 0.0;
      md->sprob[i] = 0.0;
      md->som[i] = 0.0;  
      var[i] = 0.0;
    }

    md->neval += ssamp-start;
    for (j = start; j < ssamp; j++) {
      prob = compute_prob_mixed(theta, obs, sampling, md, j, i, max_size,
				utility, 1, 0);
      /* add to the simulated probabilty */
      md->sprob[i] += (prob-md->sprob[i]) / (j+1);
//      fprintf(f, "%lg ", prob);
      /* update the second order moment */
      md->som[i] += (prob*prob-md->som[i]) / (j+1);
    }

    /* we check the value of prob to avoid failure of the
       logarithm operation */
    md->sprob[i] = nt_max(md->sprob[i], SAFE_VALUE);

    /* compute the variance */
    var[i] = sll_variance_draw(i, ssamp, md);
    *variance += var[i];

    /* gradient and likelihood contributions */
    for (k = 0; k < md->npar; k++) {
      if (md->bias)
	grad[k] += (md->sgrad[i][k]/md->sprob[i] - grad[k]
		    + (md->somgrad[i][k]/(2*md->sprob[i]*md->sprob[i])
		       - md->sgrad[i][k]*(var[i]+1)/md->sprob[i])/ssamp)/(i+1);
      else
	grad[k] += (md->sgrad[i][k]/md->sprob[i] - grad[k])/(i+1);
    }
    ll += (log(md->sprob[i])-ll)/(i+1);

//    fprintf(f, "\n");
  }

  *variance /= ssamp-1;
  *err = md->err_quantile/md->np*sqrt(*variance);
  *variance *= ssamp;

  if (md->bias)
    ll -= simulation_bias(md->np, md->err_quantile, *err);

  free(var);
  free(utility);

//  fclose(f);

  return -ll;
}

/**
 * Compute the opposite of the log-likelihood function and evaluate
 * its accuracy.
 *
 * @param theta is the vector of model parameters.
 * @param choices is a constant integer array describing the
 *   individual choices. More precisely, choices[i] is the index of
 *   the alternative chosen by individual i.
 * @param md is the data structure about the problem.
 * @param err is a pointer to a double whose content is undefined on
 *   input, and corresponds to the estimated accuracy on output, at
 *   the level md->level.
 */
double simulated_likelihood_accuracy(double *theta,
				     double ***obs, const int *choices,
				     double **sampling, const int *btype,
				     int ssamp, int max_size,
				     dcl_data *md, double *err,
				     double *variance, int start)
{
  int i, j, k;
  double ll;
  double *utility;
  double iprob;

  utility = _mallocs(md->nalt, double);

  ll = 0.0;
  *variance = 0.0;

  for (i = 0; i < md->np; i++) {

    for (k = 0; k < md->nalt; k++) utility[k] = 0.0;

    if (start == 0) {
      md->sprob[i] = 0.0;
      md->som[i] = 0.0;
    }

    for (j = start; j < ssamp; j++) {
      iprob = compute_prob_mixed(theta, obs, sampling, md, j, i,
				 max_size, utility, 0, 0);

      /* add to the simulated probabilty */
      md->sprob[i] += (iprob-md->sprob[i]) / (j+1);
		/* update the second order moment */
		md->som[i] += (iprob*iprob-md->som[i]) / (j+1);
    }

    /* we check the value of prob to avoid failure of the
       logarithm operation */
    md->sprob[i] = nt_max(md->sprob[i], SAFE_VALUE);

    /* compute the variance */
    *variance += sll_variance_draw(i, ssamp, md);

    /* likelihood contribution */
    ll += (log(md->sprob[i])-ll)/(i+1);

    /* */
    for (k = 0; k < md->nalt; k++) utility[k] = 0.0;
    md->sprob[i] = 0.0;
    md->som[i] = 0.0;
    for (j = start; j < ssamp; j++) {
      iprob = compute_prob_mixed(theta, obs, sampling, md, j, i,
				 max_size, utility, 0, 1);
      md->sprob[i] += (iprob-md->sprob[i]) / (j+1);
      md->som[i] += (iprob*iprob-md->som[i]) / (j+1);
    }
    md->sprob[i] = nt_max(md->sprob[i], SAFE_VALUE);
    
    *variance += sll_variance_draw(i, ssamp, md);
  }
	
  *variance *= ((double)ssamp/(double)(ssamp-1));
  *err = md->err_quantile/md->np*sqrt(*variance/ssamp);

  free(utility);
  if (md->bias) ll -= simulation_bias(md->np, md->err_quantile, *err);

  return -ll;
}

/**
 * WRAPPERS
 */

/**
 * @param theta is the vector of model parameters.
 * @param n is the size of theta, i.e. the number of parameters.
 * @param data is a pointer to the data structure, of type dcl_data,
 *   describing the problem under investigation.
 */
double wrapper(double *theta, int n, void *data)
{
  int shift = 0;
  dcl_data *md = (dcl_data *)data;

  return likelihood(theta, md->nobs,
		    md->obs, md->choices, md->availability,
		    md->nvar, md, shift);
}

/**
 * @param theta is the vector of model parameters.
 * @param n is the size of theta, i.e. the number of parameters.
 * @param data is a pointer to the data structure, of type dcl_data,
 *   describing the problem under investigation.
 */
double wrapper_grad(double *theta, int n, void *data, double *grad)
{
  int shift = 0;
  dcl_data *md = (dcl_data *)data;

  double ll, mean;

  ll = lg(theta, md->nobs, md->obs, md->choices, md->availability,
	  md->nvar, grad, md, shift);
  mean = ll;

  return mean;
}

/**
 * Compute the mixed-logit function.
 *
 * @param theta is the vector of model parameters.
 * @param n is the size of theta, i.e. the number of parameters.
 * @param data is a pointer to the data structure, of type dcl_data,
 *   describing the problem under investigation.
 */
double wrapper_simulated(double *theta, int n, void *data)
{
  dcl_data *md = (dcl_data *)data;

  return simulated_likelihood(theta, md->obs, md->choices,
			      md->sampling, md->nvar,
			      md->size_sampling,
			      md->max_size, md);
}

/**
 * Compute the mixed-logit function and its gradient.
 *
 * @param theta is the vector of model parameters.
 * @param n is the size of theta, i.e. the number of parameters.
 * @param data is a pointer to the data structure, of type dcl_data,
 *   describing the problem under investigation.
 */
double wrapper_simulated_grad(double *theta, int n, void *data,
			      double *grad)
{
  dcl_data *md = (dcl_data *)data;

  return slg(theta, md->obs, md->choices, md->sampling, md->btype,
	     md->size_sampling, md->max_size, md, NULL, grad);
}

/**
 * @param theta is the vector of model parameters.
 * @param n is the size of theta, i.e. the number of parameters.
 * @param data is a pointer to the data structure, of type dcl_data,
 *   describing the problem under investigation.
 *
 * @return the opposite of the log-likelihood value at theta.
 */
double wrapper_simulated_accuracy(double *theta, int n, void *data,
				  double *epsilon, double *variance,
				  int size)
{
  dcl_data *md = (dcl_data *)data;

  return simulated_likelihood_accuracy(theta, md->obs, md->choices,
				       md->sampling, md->btype,
				       md->size_sampling,
				       md->max_size, md,
				       epsilon, variance, size);
}

/**
 * @param theta is the vector of model parameters.
 * @param n is the size of theta, i.e. the number of parameters.
 * @param data is a pointer to the data structure, of type dcl_data,
 *   describing the problem under investigation.
 */
double wrapper_simulated_accuracy_grad_qmc(double *theta, int n, void *data,
					   double *epsilon, double *variance,
					   double *grad, int start)
{
  double sll, *slls;
  dcl_data *md = (dcl_data *)data;
  double mean, var;

  int shift;

  int i;

  slls = _mallocs(md->qmcrep, double);
  shift = md->max_size*md->np;


  sll = slls[0] =  slg(theta, md->obs, md->choices, md->sampling,
		       md->btype, md->size_sampling, md->max_size, md,
		       epsilon, grad);

  for (i = 1; i < md->qmcrep; i++) {
    slls[i] = simulated_likelihood(theta, md->obs, md->choices,
				   md->sampling+i*shift, md->nvar,
				   md->size_sampling, md->max_size, md);
    /* compute the variance over the various estimation */
    st_mean_variance(slls, 1, 0, md->qmcrep, &mean, &var);
  }
  *epsilon = md->err_quantile*sqrt(var);

  free(slls);

  return sll;
}

/**
 * @param theta is the vector of model parameters.
 * @param n is the size of theta, i.e. the number of parameters.
 * @param data is a pointer to the data structure, of type dcl_data,
 *   describing the problem under investigation.
 */
double wrapper_simulated_accuracy_grad(double *theta, int n, void *data,
				       double *epsilon, double *variance,
				       double *grad, int start)
{
  dcl_data *md = (dcl_data *)data;
  double ll;

  switch (md->method) {
  case MOCA_MONTE_CARLO:
    ll = sldag(theta, md->obs, md->choices, md->sampling,
		 md->size_sampling, md->max_size, md,
		 epsilon, variance, grad, start);
    return ll;
  default:
    return wrapper_simulated_accuracy_grad_qmc(theta, n, data, epsilon,
					       variance, grad, start);
  }
}
