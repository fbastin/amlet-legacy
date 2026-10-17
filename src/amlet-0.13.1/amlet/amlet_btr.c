/**
 * amlet_btr.c
 *
 * AMLET (Another Mixed Logit Estimation Tool)
 * Copyright (C) 2002-2008 Fabian Bastin <fbas@math.fundp.ac.be>
 *
 * Licensed under the MIT License
 */

#include "amlet_btr.h"

#include <math.h>
#include <float.h>
#include <limits.h>
#include <string.h>

#include <oratio/oratio.h>

#include <amlet/amlet.h>
#include "amlet_config.h"
#include "amlet_stopping.h"
#include "bhhh.h"

#define MAXNLAST    20

#define RATIO_THRESHOLD 0.2

void dcl_data_clone_btrda(dcl_data *dest, dcl_data *src)
{
  /* TODO MODIFY FOR A BETTER MEMORY GESTION */

  dest->twopowers = twopowers;

  dest->obs = src->obs;
  dest->rep = src->rep;
  dest->availability = src->availability;

  dest->choices = src->choices;
  dest->sampling = src->sampling;
  dest->size_sampling = src->size_sampling;
  dest->max_size = src->max_size;

  dest->n_random_coeff = src->n_random_coeff;

  dest->sprob = _mallocs(src->np, double);
  dest->som = _mallocs(src->np, double);
  dest->sgrad = nt_matrix_new(src->np, src->npar);
  if (src->bias) dest->somgrad = nt_matrix_new(src->np, src->npar);
  dest->bias = src->bias;

  dest->fixed = _mallocs(src->nfixed, double);
  dest->nfixed = src->nfixed;
  memcpy(dest->fixed, src->fixed, src->nfixed*sizeof(double));

  dest->nobs = src->nobs;
  dest->np = src->np;
  dest->nalt = src->nalt;
  dest->nvar = src->nvar;
  dest->npar = src->npar;

  dest->n = src->n;
  dest->btype = src->btype;
  dest->err_quantile = src->err_quantile;

  dest->np_list = src->np_list;

  dest->method = src->method;
  dest->qmcrep = src->qmcrep;

  dest->neval = src->neval;
}

/**
 * Give the candidate sampling size, following the known
 * informations, needed to achieve precision given by precision.
 *
 * @param n is the number of individuals
 *        precision: > 0;
 *        ssamp: current sampling size
 */
int candidate_size(int n, double var,
			  double quantile, double precision,
			  int ssamp, int min_size)
{
  double t;
  int size;

  t = quantile/(n*precision);
  t = t*t*var;
  size = (t > INT_MAX) ? INT_MAX : nt_max(min_size, ceil(t));

  return size;
}

void stupid_uss(double ratio, double rsize, dcl_data *md,
		       int suggested_size)
{
  /* only for debugging purposes */
  md->size_sampling = md->max_size;
}

int update_sampling_size(double ratio, double rsize,
			  int suggested_size, int min_size, int max_size)
{
  int size_sampling;
  const double a1 = 0.5, a2 = 1.0;

  /*  if (verbose)
    printf("Update sampling size: ratio %g, rsize %g, suggested size %d\n",
    ratio, rsize, suggested_size); */

  if (ratio >= 1.0)
    size_sampling = nt_min(a1*max_size, a2*suggested_size);
  else if (ratio >= rsize)
    /* the ratio between the observed decrease and the precision
       is greater the the ratio between the used size and the
       suggested size */
    size_sampling = nt_min(a1*max_size, ratio*suggested_size);
  else if (ratio >= RATIO_THRESHOLD)
    size_sampling = nt_min(a1*max_size, suggested_size);
  else
    size_sampling = max_size;

  size_sampling = nt_max(min_size, size_sampling);

  /*  if (verbose)
      printf("End update sampling size: %d\n", size_sampling); */

  return size_sampling;
}


int update_sampling_size_rtr(double ratio, double rsize, int current,
			     int suggested_size, int min_size, int max_size)
{
  int size_sampling;
  const double a1 = 0.5, a2 = 1.0;

  /*  if (verbose)
    printf("Update sampling size: ratio %g, rsize %g, suggested size %d\n",
    ratio, rsize, suggested_size); */

  if (ratio >= 1.0)
    size_sampling = nt_min(a1*max_size, a2*suggested_size);
  else if (ratio >= rsize)
    /* the ratio between the observed decrease and the precision
       is greater the the ratio between the used size and the
       suggested size */
    size_sampling = nt_min(a1*max_size, ratio*suggested_size);
  else if (ratio >= RATIO_THRESHOLD)
    size_sampling = nt_min(a1*max_size, suggested_size);
  else
    size_sampling = max_size;

  size_sampling = nt_min(size_sampling, 2*current);

  size_sampling = nt_max(min_size, size_sampling);

  /*  if (verbose)
      printf("End update sampling size: %d\n", size_sampling); */

  return size_sampling;
}

/**
 * BTR with dynamic accuracy
 *
 * Ref.:
 * F. Bastin, C. Cirillo, and Ph. L. Toint, An adaptive Monte Carlo
 * algorithm for computing mixed logit estimators, Computational
 * Management Science, 3(1):55-79, 2006.
 *
 * @param f is the objective function.
 * @param fg is the objective and gradient.
 * @param tol is the user-preset tolerance.
 * @param btr is btr object (must be initialized).
 *
 * @return the optimum value
 */
double model_amlet_btr_da(double (*f)(double *, int, void *, double *,
				      double *, int),
			  double (*fg)(double *, int, void *, double *,
				       double *, double *, int),
			  void *data, double tol, BTR *btr,
			  double **hessian,
			  const enum CBLAS_ORDER Order,
			  const enum CBLAS_UPLO Uplo,
			  double *work, double *epsilon)
{
  int n;  /* number of parameters */

  double new_value, rho = 0;
  double *y = NULL, *newgrad = NULL;
  double *x_cand;
  double *oldvals;
  int    *last;

  int nan = 0;

  double model_decrease = 0.0;  /* decrease in the model */
  double fct_decrease = 0.0;    /* decrease in the function */

  double epsilonplus;
  double variance, variance2;

  int nlast;

  double rfct, rmodel, rsize, tmp;

  dcl_data *md = (dcl_data *)data;  /* current iterate */
  dcl_data *md2;                      /* next iterate */
  dcl_data *mdtmp;

  int i;
  int start_ssize, suggested_size = 0;
  int tmp_ssize, min_size, newsize;

  int success = 0;   /* counter of successful iterates */

  double dval;
  int dlast;

  /* Initialization */
  n = btr->n;
  btr->k = 0; /* set the iteration counter */

  nlast = 0;

  /* initial sampling size */
  min_size = MINSAMPLING;
  md->size_sampling = nt_max(min_size, rint(md->max_size/10.0));
  md->neval = 0;
  /* md->size_sampling = md->max_size; */

  md2 = _malloc(dcl_data);
  dcl_data_clone_btrda(md2, md);

  /* Memory allocation */
  x_cand = _mallocs(n, double);
  oldvals = _mallocs(md->max_size, double);
  oldvals--;  /* index 1 now corresponds to the first element */
  last = _mallocs(md->max_size, int);
  last--;
  if (btr->d2approx || fg) {
    y = _mallocs(n, double);
    newgrad = _mallocs(n, double);
  }
  else
    y = _mallocs(n, double);

  start_ssize = 0;
  if (fg) {
    btr->value = fg(btr->x, n, md, epsilon, &variance, btr->grad, 0);
  }
  else {
    btr->value = f(btr->x, n, md, epsilon, &variance, 0);
    (btr->d1)(btr->x, n, md, btr->grad);
  }
  oldvals[md->size_sampling] = btr->value;
  last[md->size_sampling] = 0;

  for (i = 0; i < n; i++) btr->step[i] = DBL_MAX;
  rmodel = rfct = DBL_MAX;
  for (i = 1; i <= md->max_size; i++) {
    oldvals[i] = DBL_MAX;
    last[i] = -1;
  }

  /* Use only for debugging purposes */
  /*  nt_vect_print(stdout, "g", btr->grad, n);
  nt_derive_grad_fd(wrapper_simulated,
		    btr->x, btr->grad, work, 0.0, work+n, n, NULL, work+2*n,
		    data);
		    nt_vect_print(stdout, "g", btr->grad, n); */

  /* Review stopping test? */
  while (!(btr->info = amlet_btrda_stopping_test(btr, btr->step, nlast,
						 md, *epsilon, rfct, rho))
	 && !btr->info) {

    /*    nt_vect_print(stdout, "g", btr->grad, n);
    nt_derive_grad_fd(wrapper_simulated,
		      btr->x, btr->grad, work, 0.0, work+n, n, NULL, work+2*n,
		      data);
		      nt_vect_print(stdout, "g", btr->grad, n); */
    if (btr->verbose) printf("Iteration %d\n", btr->k);

    if (btr->printer) btr->printer(btr->output, btr, md);   
    tmp_ssize = md->size_sampling;

    /* step calculation */
    btr->step_function(btr, CblasRowMajor, CblasUpper, hessian, work,
		       data);

    model_decrease = model(CblasRowMajor, CblasUpper, btr->step,
			   btr->grad, *hessian, n, data);

    /* compute the candidate sampling size */
    suggested_size = candidate_size(md->np, variance,
				    md->err_quantile, model_decrease,
				    md->size_sampling, min_size);

    rmodel = (*epsilon > 0.0) ? model_decrease/(*epsilon) : DBL_MAX;
    rsize = (double)md->size_sampling / nt_min(md->max_size, suggested_size);

    md2->size_sampling = update_sampling_size(rmodel, rsize,
					      suggested_size, min_size,
					      md2->max_size);

    /* stupid_uss(rmodel, rsize, md, suggested_size); */

    cblas_dcopy(n, btr->x, 1, x_cand, 1);
    cblas_daxpy(n, 1.0, btr->step, 1, x_cand, 1);

    if (fg)
      new_value = fg(x_cand, n, md2, &epsilonplus, &variance2,
		     newgrad, 0);
    else
      new_value = f(x_cand, n, md2, &epsilonplus, &variance2, 0);

#ifdef NAN
    if (isnan(new_value) || isinf(new_value)) {
      nan = 1;
      if (btr->verbose)
	fprintf(stderr,
		"Incorrect function value; reduce the trust-region.\n");
      nt_vect_print(stdout, "xc", x_cand, n);
      printf("new_value = %g\n", new_value);
      rho = 0;
      /* Not very good, but shoult not occur in practice. */
      md2->size_sampling = 36;
    }
    else {
#endif
      /* acceptance of the trial point */
      fct_decrease = btr->value-new_value;
      rho = fct_decrease/model_decrease;

      if (rho < btr->eta1 && tmp_ssize != md2->size_sampling) {

	/* Perhaps we could use the following code, in conjunction
	   with bias estimation */
	/* if (fabs(fct_decrease) <= nt_max(*epsilon, epsilonplus))
	  newsize = ceil(0.5*(md->size_sampling+md2->size_sampling));
	  else */
	if (tmp_ssize < md2->size_sampling) {
	  /* the sampling size has been increased */
	  start_ssize = tmp_ssize;
	  md->size_sampling = md2->size_sampling;
	  if (fg)
	    btr->value = fg(btr->x, n, md, epsilon, &variance,
			    btr->grad, tmp_ssize);
	  else {
	    btr->value = f(btr->x, n, md, epsilon, &variance, tmp_ssize);
	    /* could be optimized */
	    (btr->d1)(btr->x, n, md, btr->grad);
	  }

	  /* Updating the model_decrease, as initially proposed, could
	     produce problems, since nothing ensures that the computed
	     step is still valid for the model using the updated size.
	     The following line is still present to remember this fact.*/
	  /* model_decrease = model(CblasRowMajor, CblasUpper, btr->step,
	     btr->grad, hessian, n, md); */

	  fct_decrease = btr->value-new_value;
	  rho = fct_decrease/model_decrease;

	  if (btr->printer)
	    printf("Corrected sampling size: %d\n", md2->size_sampling);
	}
	else {
	  /* the sampling size has been decreased */
	  i = 0;
	  /* Compute the sample size correspondinf to the bias */
	  tmp = 0.5/(md->np*model_decrease)*variance;
	  newsize = (tmp > INT_MAX) ?
	    INT_MAX : nt_max(md2->size_sampling, ceil(tmp));
	  if (newsize == md2->size_sampling)
	    newsize = md->size_sampling;
	  else {
	    newsize = nt_min(md->size_sampling, newsize);
	    //	    printf("RESULAT %d\n", newsize);
	  }
	  //	  newsize = md->size_sampling;

	  do {
	    start_ssize = md2->size_sampling;
	    md2->size_sampling = newsize;
	    if (fg)
	      new_value = fg(x_cand, n, md2, &epsilonplus, &variance2,
			     newgrad, start_ssize);
	    else
	      new_value = f(x_cand, n, md2,
			    &epsilonplus, &variance2, start_ssize);

	    fct_decrease = btr->value-new_value;
	    rho = fct_decrease/model_decrease;
	    if (rho >= btr->eta1 || newsize == md->size_sampling)
	      i = 1;
	    else
	      newsize = md->size_sampling;
	  } while (!i);
	}

	rmodel = (*epsilon > 0.0) ? model_decrease/(*epsilon) : DBL_MAX;
      }
      rfct = (*epsilon > 0.0) ? fct_decrease/(*epsilon) : DBL_MAX;
    
      /* If Hessian approximation is used we compute
	 the gradients difference. */
      if (btr->d2approx) {
	if (!fg)
	  (btr->d1)(x_cand, n, md2, newgrad);
	cblas_dcopy(n, newgrad, 1, y, 1);
	cblas_daxpy(n, -1.0, btr->grad, 1, y, 1);
	if (!nan)
	  ((deriv2_approx)(btr->d2approx))(CblasRowMajor, CblasUpper, n,
					   *hessian, n, y, 1,
					   btr->step, 1, work);
	nan = 0.0;
      }
#ifdef NAN
    }
#endif

#ifdef DEBUG
    printf("rfct: %g, rmodel: %g epsilon: %g\nmodel_decrease: %g fct_decrease: %g rho: %g\n\n ", rfct, rmodel, *epsilon, model_decrease, fct_decrease, rho);
#endif

    if (rho >= btr->eta1) ++success;

    if (tmp_ssize != md2->size_sampling) {
      if (rho >= btr->eta1) {
	dlast = success-last[md2->size_sampling];
	dval = oldvals[md2->size_sampling]-new_value;
	if (dval < 0.5*RATIO_THRESHOLD*dlast*epsilonplus) {
	  if (tmp_ssize < md2->size_sampling) {
	    /* we have increased the sampling size */
	    min_size = nt_min(ceil(0.5*(md2->size_sampling+tmp_ssize)),
			    md2->max_size);
	    /*	    min_size = nt_min(tmp_ssize+
		    ceil(M_SQRT1_2*(md2->size_sampling-tmp_ssize)),
		    md2->max_size); */
	  }
	  else
	    /* we have decreased the sampling size */
	    min_size = nt_min(md2->size_sampling+1, md2->max_size);
	}

	/* record the value */
	oldvals[md2->size_sampling] = new_value;
	last[md2->size_sampling] = success;
      }
      else if (tmp_ssize != md->size_sampling) {
	dlast = success-last[md->size_sampling];
	dval = oldvals[md->size_sampling]-btr->value;
	if (dval < 0.5*RATIO_THRESHOLD*dlast*(*epsilon))
	  min_size = nt_min(ceil(0.5*(md2->size_sampling+tmp_ssize)),
			  md->max_size);
	/* record the value */
	oldvals[md->size_sampling] = btr->value;
	last[md->size_sampling] = success;
      }

      if (btr->verbose) printf("Minimum sample size: %d\n", min_size);
    }

    /* Update the trust region radius */
    btr->update_radius(btr, rho);

    /* Update the iterate */
    if (nlp_btr_accept_candidate(btr, rho)) {
      /* Set the new iterate to the candidate point. */
      cblas_dcopy(btr->n, x_cand, 1, btr->x, 1);
      btr->value = new_value;

      *epsilon = epsilonplus;
      variance = variance2;
      mdtmp = md;
      md = md2;   /* needed to have the correct information on partial
		     gradients and variance */
      md2 = mdtmp;

      if (rfct >= 0.01)
	nlast = 0;
      else if (rho >= btr->eta2) nlast++;

      if (fg || btr->d2approx)
	/* we already know the gradient */
	cblas_dcopy(n, newgrad, 1, btr->grad, 1);
      else {
	cblas_dcopy(n, btr->grad, 1, y, 1);
	(btr->d1)(btr->x, n, md, btr->grad);
      }

      if (!btr->d2approx) {
	cblas_daxpy(n, -1.0, btr->grad, 1, y, 1);
        cblas_dscal(n, -1.0, y, 1);
      }
    }

    if ((++(btr->k)) == btr->max_iter)
      btr->info = NLP_STOPPING_ITERATIONS;
  }

  /* print results */
  btr_print_results(btr->output, btr, md);

  /* memory deallocation */

  free(x_cand);
  free(++oldvals);
  free(++last);
  if (y) free(y);
  free(newgrad);

  printf("NEVAL: %d\n", md->neval);
  /* BUG!?!?!? */
  /*  if (data == md2) {
    mdtmp = md;
    md = md2;
    md2 = mdtmp;
  }
  free(md2->sprob);
  free(md2->som);
  nt_matrix_free(md2->sgrad, md2->np, md2->npar);
  free(md2);*/

  return btr->value;
}

/**
 * BTR with dynamic accuracy
 *
 * Ref.:
 * F. Bastin, C. Cirillo, and Ph. L. Toint, An adaptive Monte Carlo
 * algorithm for computing mixed logit estimators, Computational
 * Management Science, 3(1):55-79, 2006.
 *
 * @param f is the objective function.
 * @param fg is the objective and gradient.
 * @param acc is the accuracy evaluation function.
 * @param tol is the user-preset tolerance.
 * @param btr is btr object (must be initialized).
 *
 * @return the optimum value
 */
double model_amlet_rtr_da(double (*f)(double *, int, void *, double *,
				      double *, int),
			  double (*fg)(double *, int, void *, double *,
				       double *, double *, int),
			  void *data, double tol, BTR *btr,
			  double **hessian,
			  const enum CBLAS_ORDER Order,
			  const enum CBLAS_UPLO Uplo,
			  double *work, double *epsilon)
{
  int n;  /* number of parameters */

  double new_value[32], rho = 0;
  double *y = NULL;
  double *x_cand;
  double *oldvals;
  int    *last;

  int nan = 0;

  double model_decrease = 0.0;  /* decrease in the model */
  double fct_decrease = 0.0;    /* decrease in the function */

  double epsilonplus[32];
  double variance, variance2[32];

  int nlast;

  double rfct, rmodel, rsize;

  dcl_data *md = (dcl_data *)data;  /* current iterate */
  int new_size;

  int i;
  int invert;
  int start_ssize, suggested_size = 0;
  int tmp_ssize, min_size;

  int success = 0;   /* counter of successful iterates */

  FILE *log;

  int pstart, restart;

  log = fopen("sizes.txt", "w");

  /* Initialization */
  n = btr->n;
  btr->k = 0; /* set the iteration counter */
  btr->retro = 1;
  md->neval = 0;

  nlast = 0;

  /* initial sampling size */
  min_size = MINSAMPLING;
  md->size_sampling = min_size;//nt_max(min_size, rint(md->max_size/10.0));
  /* md->size_sampling = md->max_size; */

  /* Memory allocation */
  x_cand = _mallocs(n, double);
  oldvals = _mallocs(md->max_size, double);
  oldvals--;  /* index 1 now corresponds to the first element */
  last = _mallocs(md->max_size, int);
  last--;
  y = _mallocs(n, double);

  start_ssize = 0;
  if (fg) {
    btr->value = fg(btr->x, n, md, epsilon, &variance, btr->grad, 0);
  }
  else {
    btr->value = f(btr->x, n, md, epsilon, &variance, 0);
    (btr->d1)(btr->x, n, md, btr->grad);
  }
  oldvals[md->size_sampling] = btr->value;
  last[md->size_sampling] = 0;

  for (i = 0; i < n; i++) btr->step[i] = DBL_MAX;
  rmodel = rfct = DBL_MAX;
  for (i = 1; i <= md->max_size; i++) {
    oldvals[i] = DBL_MAX;
    last[i] = -1;
  }

  /* Review stopping test? */
  while (!(btr->info = amlet_btrda_stopping_test(btr, btr->step, nlast,
						 md, *epsilon, rfct, rho))
	 && !btr->info) {

    invert = 0;
    if (btr->verbose) printf("Iteration %d\n", btr->k);

    if (btr->printer) btr->printer(btr->output, btr, md);   
    tmp_ssize = md->size_sampling;

    /* step calculation */
    btr->step_function(btr, CblasRowMajor, CblasUpper, hessian, work,
		       data);

    model_decrease =
      btr_quadratic_model(CblasRowMajor, CblasUpper, n, btr->step,
			  btr->grad, *hessian, data, work);
    cblas_dcopy(n, btr->x, 1, x_cand, 1);
    cblas_daxpy(n, 1.0, btr->step, 1, x_cand, 1);

    /* Compute VALUES */
    /* this should be iterative to store values */
    restart = 0;
    pstart = 4;

    md->size_sampling = md->twopowers[pstart];
    do {
      md->size_sampling = nt_min(tmp_ssize, md->twopowers[pstart]);
      if (fg)
	new_value[pstart] = fg(x_cand, n, md, epsilonplus+pstart,
			       variance2+pstart, md->newgrad[pstart], restart);
      else
	new_value[pstart] = f(x_cand, n, md, epsilonplus+pstart,
			      variance2+pstart, restart);
      restart = md->twopowers[pstart++];
    } while (md->size_sampling < tmp_ssize);
    pstart--;

    /* acceptance of the trial point */
    fct_decrease = btr->value-new_value[pstart];
    rho = fct_decrease/model_decrease;
    //    new_value[0] = new_value[pstart-1];

    /* If Hessian approximation is used we compute
       the gradients difference. */
    if (btr->d2approx) {
      if (!fg)
	(btr->d1)(x_cand, n, md, md->newgrad[pstart]);
      cblas_dcopy(n, md->newgrad[pstart], 1, y, 1);
      cblas_daxpy(n, -1.0, btr->grad, 1, y, 1);
      if (!nan)
	((deriv2_approx)(btr->d2approx))(CblasRowMajor, CblasUpper, n,
					 *hessian, n, y, 1,
					 btr->step, 1, work);
      nan = 0.0;
    }

    rfct = (*epsilon > 0.0) ? fct_decrease/(*epsilon) : DBL_MAX;
    rmodel = (*epsilon > 0.0) ? model_decrease/(*epsilon) : DBL_MAX;

    model_decrease =
      btr_quadratic_retro_model(CblasRowMajor, CblasUpper, n, btr->step,
				md->newgrad[pstart], *hessian, data, work);

    /* compute the candidate sampling size */
    suggested_size = candidate_size(md->np, variance2[pstart],
				    md->err_quantile, fct_decrease,
				    md->size_sampling, min_size);
      
    rsize = (double)md->size_sampling / nt_min(md->max_size, suggested_size);

    new_size = update_sampling_size_rtr(rmodel, rsize, tmp_ssize, suggested_size, min_size, md->max_size);

    // tester la norme du gradient!!!

    printf("%d %g %d %d\n", btr->k, rho, tmp_ssize, new_size);
    fprintf(log, "%d %d\n", btr->k, tmp_ssize);

    if (rho >= btr->eta1 && tmp_ssize != new_size) {
      /* The candidate iterate will be accepted. */
      /* THIS CAN BE IMPROVED! WE COMPUTE TOO MUCH INFORMATION WHEN THE SAMPLE SIZE IS DECREASED! */
      start_ssize = tmp_ssize;
      if (tmp_ssize <= new_size) {
	md->size_sampling = new_size;
	//      printf("START: %d\n", start_ssize);
	if (fg)
	  new_value[pstart] = fg(x_cand, n, md, epsilonplus+pstart,
				 variance2+pstart, work, start_ssize);
	else
	  new_value[pstart] = f(x_cand, n, md,
			   epsilonplus+pstart, variance2+pstart, start_ssize);
	//      if (test_value < btr->value || new_size > tmp_ssize) {
	cblas_dcopy(n, work, 1, md->newgrad[pstart], 1);
      }
      else {
	/* get the value back */
      	pstart = nt_two_power_ceil(new_size);
	md->size_sampling = nt_min(twopowers[pstart], md->max_size);
	/* if md->size_sampling is greater than tmp_ssize, we have not
	   computed the value... we cannot increase the sample size! */
	md->size_sampling = nt_min(md->size_sampling, tmp_ssize);

	/* Only the following line for the monotonic version. */
	//	md->size_sampling = tmp_ssize;

      }
      invert = 1;
    }
 
#ifdef DEBUG
    printf("rfct: %g, rmodel: %g epsilon: %g\nmodel_decrease: %g fct_decrease: %g rho: %g\n\n ", rfct, rmodel, *epsilon, model_decrease, fct_decrease, rho);
#endif

    // AJOUTER UN TEST DE TAILLE MINIMALE

    if (rho >= btr->eta1) ++success;

    /* Update the iterate */
    if (nlp_btr_accept_candidate(btr, rho)) {
      /* Set the new iterate to the candidate point. */
      cblas_dcopy(btr->n, x_cand, 1, btr->x, 1);
      btr->value = new_value[pstart];

      *epsilon = epsilonplus[pstart];
      variance = variance2[pstart];

      if (rfct >= 0.01)
	nlast = 0;
      else if (rho >= btr->eta2) nlast++;

      if (fg || btr->d2approx)
	/* we already know the gradient */
	cblas_dcopy(n, md->newgrad[pstart], 1, btr->grad, 1);
      else {
	cblas_dcopy(n, btr->grad, 1, y, 1);
	(btr->d1)(btr->x, n, md, btr->grad);
      }

      if (!btr->d2approx) {
	cblas_daxpy(n, -1.0, btr->grad, 1, y, 1);
        cblas_dscal(n, -1.0, y, 1);
      }
    }

    if (btr->retro && (rho >= btr->eta1)) {
      rho = fct_decrease /
	btr_quadratic_retro_model(Order, Uplo, n, btr->step,
				  btr->grad, *hessian, data, work);
    }

    /* Update the trust region radius */
    btr->update_radius(btr, rho);

    if ((++(btr->k)) == btr->max_iter)
      btr->info = NLP_STOPPING_ITERATIONS;
  }

  /* print results */
  btr_print_results(btr->output, btr, md);

  /* memory deallocation */

  free(x_cand);
  free(++oldvals);
  free(++last);
  if (y) free(y);

  printf("NEVAL: %d\n", md->neval);
  fclose(log);

  return btr->value;
}
