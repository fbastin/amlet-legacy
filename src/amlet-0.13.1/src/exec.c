/**
 * amlet.c
 *
 * AMLET (Another Mixed Logit Estimation Tool)
 * Copyright (C) 2002-2012 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include "amlet/amlet.h"
#include "info.h"
#include "amlet/amlet_config.h"

#include <getopt.h>

#include <oratio/oratio.h>
#include <amlet/dcl.h>
#include <ophelia/solutions.h>

#include <oratio/ctimer.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "alternative.h"
#include "amlet_analysis.h"
#include "amlet_data.h"
#include "amlet_stopping.h"
#include "bhhh.h"
#include "bootstrap.h"
#include "amlet_driver.h"
#include "amlet/external_strategy.h"
#include "io.h"
#include "log.h"
#include "observations.h"
#include "amlet/optimization.h"

#include "amlet_btr.h"

#include <config.h>

/* TODO: check how to properly free the memory in case of failure do to detected error */

/**
 * Set the simulation parameters.
 *
 * @param m is the model to estimate.
 * @param par on output, the value of parameters used for the
 *   generation process
 */
int set_parameters(AMLET_Model *m,
		   amlet_config *ac, AMLET_Driver *dr, int nvar,
		   int *npar, int *n_random_coeff, int **ff, int *nfixed,
		   int **btype, double **x, double **fixed)
{
  int i, j, k, l, t;
  void *np_object;   /* nonparametric object */
  DS_List *np_list = data.np_list;
  void *obj;
  int s = 0;
  //  OP_Problem *p = amlet_get_problem(m);

  np_list = ds_list_new();

  *n_random_coeff = 0;
  *nfixed = 0;

  if (dr->array[BTYPE].n != nvar)
    amlet_error(ac, mclog->f, 1, AMLET_ERROR_PARSING,
		"incorrect number of parameters in btype declaration");

  /* We set the type vector, following the informations of the driver. */
  *btype = _mallocs(nvar, int);
  driver_affect_integer(*btype, dr->array+BTYPE);

  /* We make a first screening to determine the parameters types, and
     their numbers. */
  for (i = 0; i < nvar; i++) {
    np_object = NULL;
    t = (*btype)[i];
    if (t != CONSTANT && t != FIXED) (*n_random_coeff)++;
    /* if nonparametric... we have to register the right object */
    if (amlet_is_nonparametric(t)) {
      np_object = amlet_nonparametric_set(t);
      np_list = ds_list_prepend(np_list, np_object);
    }
    if (t == FIXED)
      (*nfixed)++;
    else
      s += moca_count_parameters(t, np_object);
  }

  np_list = ds_list_reverse(np_list);

  driver_free_array(dr, BTYPE);  /* SHOULD BE FREED, BUT USED IN MCBTRC */

  np_list = ds_list_first(np_list);
  data.np_list = np_list;

  /* We make a second screening to affect the initial parameters
     values. */
  *npar = s;
  data.beta = _mallocs(*npar, double);
  *x = _mallocs(s, double);
  *fixed = _mallocs(*nfixed, double);
  *ff = _mallocs((*npar+*nfixed), int);

  /* WE SHOULD ALSO TEST IF NFIXED IS EQUAL TO 0. */
  if (!dr->array[B].array) {
    amlet_ml_starting_point(*npar, *x, 1);
    return 0;
  }

  if (dr->array[B].n != *npar + *nfixed)
    return amlet_error(ac, mclog->f, 1, AMLET_ERROR_STARTING_POINT,
		       "incorrect number of parameters in b declaration %d; should be %d", dr->array[B].n, *npar);
  
  j = k = l = t = 0;
  for (i = 0; i < (s+*nfixed); i++) (*ff)[i] = 0;

  for (i = 0; i < nvar; i++) {
    switch ((*btype)[i]) {
    case FIXED:
      (*fixed)[t++] = (dr->array[B].array[l].type == REAL) ?
	dr->array[B].array[l].d.d : dr->array[B].array[l].d.i;
      (*ff)[l++] = 1;
      break;
    case RAN_SP_SPLINE:
      obj = np_list->data;
      np_list = np_list->next;
      for (j = 0; j < moca_count_parameters((*btype)[i], obj); j++) {
	(*x)[k++] = (dr->array[B].array[l].type == REAL) ?
	  dr->array[B].array[l].d.d : dr->array[B].array[l].d.i;
	l++;
      }
      break;
    default:
      for (j = 0; j < moca_count_parameters((*btype)[i], NULL); j++) {
	(*x)[k++] = (dr->array[B].array[l].type == REAL) ?
	  dr->array[B].array[l].d.d : dr->array[B].array[l].d.i;
	l++;
      }
    }
  }
  driver_free_array(dr, B);

  return 0;
}

/**
 * Arguments parsing
 */

/**
 * Parse command line arguments
 *
 * @return AMLET_NORMAL_EXIT code if the parsing has not encountered
 *   any error, and the corresponding code otherwise.
 */
int parse_args(AMLET_Model *m, int argc, char **argv, NTLog *log)
{
  int c, option_index;

  /* Declaration of options. */
  static struct option long_options[] = {
    {"debug", no_argument, 0, 'd'},
    {"help", no_argument, 0, 'h'},
    {"output", required_argument, 0, 'o'},
    {"verbose", no_argument, 0, 'v'},
    {"version", no_argument, 0, 'V'},
    {0, 0, 0, 0}   /* Termination of the long_options list */
  };

  while ((c = getopt_long(argc, argv, "dhlo:uvV", long_options,
			  &option_index)) != -1) {
    switch(c) {

    case 'd':
      break;
    case 'h':
      amlet_help(NULL);
      return AMLET_HELP_EXIT;
    case 'V':
      amlet_version();
      return AMLET_VERSION_EXIT;
    case 'o':
      nt_log_set_file(log, optarg);
      break;
    case 'v':
      m->verbose = 1;
      printf("Verbose mode\n");
      break;
    default:
      return AMLET_ERROR_OPTIONS;
    }
  }

  return AMLET_NORMAL_EXIT;
}

/**
 * Main program
 */
int amlet_exec(int argc, char *argv[])
{
  AMLET_Model *m;
  OP_Problem *p;
  OP_Solutions *solutions = NULL;

  FILE *src;
  FILE *f = NULL;
  AMLET_Driver dr;    /* driver for the calibration process */

  amlet_config ac;
  int info = 0; /* info identifier on AMLET behaviour; it will be set
		   to 0 if no problem has been encountered. */

  char *buff = NULL;  /* general purpose string */
  char *str = NULL;

  int n_random_coeff; /* for mixed-logit models, number of random
			 coefficients in the utility */

  int nobs;
  double ***obs = NULL, ***n_obs = NULL;
  /* obs: three-dimensional array; first dimension corresponds to the
     observation index, second to the alternative and the last to the
     variable */
  int *rep = NULL; /* number of repetion of observations per individu */
  int *n_rep = NULL; /* number of repetion of observations per individu */
  int **availability = NULL; /* alternatives availability */
  int **n_availability = NULL; /* alternatives availability */
  int *choices = NULL;    /* choices made by the individuals */
  int *n_choices = NULL;

  /* Monte Carlo / Quasi-Monte Carlo parameters */
  Random *ran = NULL;
  unsigned long seed = 0;

  int *sampling_types = NULL; /* basic types for random variables */

  double tol;  /* tolerance used in the calibration process */

  double *work = NULL;

  DS_List *tmp;

  /* external variable sample size strategy */
  AmletExternalStrategy *aes = NULL;

  /* Optimization variables */
  DclOptiMethod opti;
  double *x = NULL, *grad = NULL;

  int niter = 0;   /* number of iterations */

  int *btype = NULL;
  int npar = 0;        /* number of parameters for the calibration process */

  double bias, accuracy = -1.0; /* bias and accuracy of the likelihood for
				   mixed-logit models */
  double elapsed_time;   /* time consumed by the optimization process */

  long sec, ms;

  /* Likelihood variables */
  double ll = 0.0;   /* likelihood value */
  double ll0 = 0.0;  /* likelihood value if all parameters are set to 0 */
  double llc = 0.0;  /* likelihood value evaluated only at independant terms */

  int i, j, out;

  int np = 0;
  int *ff = NULL;
  char *drivername = NULL;

  /* PROGRAM */

  /* Credits */
  amlet_credits(stdout, "", "");

  /* Some initialization tasks */
  ctools_timer_init();

  m = amlet_model_new();
  p = amlet_get_problem(m);

  opti.b = NULL;
  ac.logname = set_logname(&buff, &dr, &ac);
  mclog = nt_log_new(buff);

  /* Arguments parsing */
  buff = NULL;
  if (argc > 1) {
    out = parse_args(m, argc, argv, mclog);
    if (out != AMLET_NORMAL_EXIT) return out;
  }
  if (optind != argc-1) {
    amlet_help("Invalid usage: you must specify the driver file");
    exit(AMLET_ERROR_OPTIONS);
  }
  drivername = argv[optind];
  if (!(src = c_fopen(drivername, "r")))
    c_error("unable to open driver file '%s'\n", drivername);

  /* We parse the data needed for the calibration */
  if (amlet_parse_driver(src, &dr)) {
    fprintf(stderr, "Please check driver file %s\n", drivername);
    info = AMLET_ERROR_PARSING;
    goto End;
  }

  driver_assign_model(&dr, m);

  data.nvar = dr.var[NVAR].i;
  niter = dr.var[NITER].i;
  if (dr.var[CENSOR].i) data.nvar--;
  if (dr.var[IDDEP].i) data.nvar--;

  data.np = dr.var[NP].i;       /* number of people */
  data.nobs = dr.var[NOBS].i;   /* total number of observations */
  data.nalt = dr.var[NALT].i;   /* number of alternatves */

  data.n = NULL;                /* nests are currently not supported. */
  data.beta = NULL;

  data.bias = dr.var[BIAS].i;

  /* Parameters array initialization */
  if (set_parameters(m, &ac, &dr, data.nvar, &npar, &n_random_coeff,
		     &ff, &(data.nfixed), &btype, &x, &(data.fixed)))
    goto End;
  /* Check with the driver */
  driver_check_parameters(&dr, &ac, &data);

  data.n_random_coeff = n_random_coeff;

  /* Random generation */
  seed = dr.var[SEED1].i;
  ran = ran_random_new(dr.var[RANGEN].i, &seed, 1);

  if (m->verbose) {
    nt_log_section(mclog, "Random generator");
    nt_log_write(mclog, "Random generator used: %s\n",
		 ran_random_generator_name(ran));
    nt_log_write(mclog, "Seed for random generator: %d\n", seed);
  }

  dcl_memory_allocate(data.nobs, data.nalt, data.nvar, dr.var[CENSOR].i,
		      &choices, &obs, &availability);

  if (data.np <= 0) dr.var[NP].i = data.np = data.nobs;
  /* Initialize the repetitions table */
  rep = _mallocs(data.np+1, int);
  rep[0] = 0;

  /* Initialize the repetitions array. We assume that if no
     repetitions file is given, each individual produces the same
     number of observations. */ 
  j = data.nobs/data.np;
  for (i = 1; i <= data.np; i++) rep[i] = j;

  /* Observations */
  if (dr.generation)
    amlet_generate_observations(&ac, ran, &dr, data.nvar,
				choices, obs, availability, rep);
  else {
    /* Read a Gauss file */
    if (m->verbose) printf("Read data.\n");
    if (!dr.obsname)
      c_error("no specified observations");

    out = import(dr.var[FORMAT].i, dr.obsname,
		 obs, availability, choices, data.np, data.nobs,
		 data.nalt, data.nvar, dr.var[IDCENSOR].i,
		 dr.var[IDDEP].i, rep);

    if (out) {
      info = out;
      goto End;
    }
  }

  /* Check the number of observations */
  info = amlet_check_repetitions(data.np, data.nobs, rep);
  if (info) goto End;

  /* Check that individuals choose available alternatives */
  if (amlet_error_choice(mclog->f, data.nobs, choices, availability)) {
    amlet_warning(stderr, 1,
		  "invalid choices; see %s for details.", AMLET_DEFAULT_LOG);
    info = AMLET_ERROR_CHOICES;
    goto End;
  }

  int nsim = 1;
  int simulation;

  nt_two_powers_array(31, twopowers);
  data.twopowers = twopowers;

  data.npar = npar;
  data.btype = btype;
  data.sgrad = NULL;
  data.somgrad = NULL;
  data.sampling = NULL;
  nobs = data.nobs;

  data.method = dr.var[SAMPLING].i;
  data.qmcrep = dr.var[QMCREP].i;
    
  /* TODO: IS err_quantile really necessary? */
  data.err_quantile = st_gaussian_icdf(1-dr.var[LEVEL_ACCURACY].d/2.0);

  work = _mallocs(nt_max(4*npar, npar+data.nfixed), double);
  amlet_model_init(m, npar);

  int nsize = (n_random_coeff) ? data.np : data.nobs;
  data.sprob = _mallocs(nsize, double);
  data.som = _mallocs(nsize, double);
  data.sgrad = nt_matrix_new(nsize, npar);
  if (data.bias) data.somgrad = nt_matrix_new(nsize, npar);
  data.newgrad = nt_matrix_new(nsize, npar);

  if (nsim > 1) solutions = op_solutions_new(nsim, npar);

  for (simulation = 0; simulation < nsim; simulation++) {
    if (m->verbose) printf("Simulation %d\n\n", simulation);

    bootstrap_resample(&data, simulation, &(data.nobs), nobs,
			obs, rep, choices, availability,
			&n_obs, &n_rep, &n_choices, &n_availability);

    nfctev = 0;  /* initialize the function evaluations counter */

    /* Optimization algorithm */
    switch(dr.var[METHOD].i) {
    case MCBTR:
    case MCRTR:
    case MCBTRC:
    case MCBTRDA:
    case MCRTRDA:
      /* Trust-region method */
      opti.b = op_btr_new();

      op_btr_init(opti.b, npar+dr.var[NNESTS].i-1, 0); 

      memcpy(opti.b->x, x, (npar+dr.var[NNESTS].i-1)*sizeof(double));
      x = opti.b->x;

      grad = opti.b->grad;

      amlet_init_btr(opti.b, &dr);

      break;
    case MCBFGS:
    case MCBFGSDA:
      /* BFGS linesearch method */

      opti.bf = _malloc(BFGS);
      nlp_bfgs_init(opti.bf, npar+dr.var[NNESTS].i-1, 0);

      memcpy(opti.bf->x, x, (npar+dr.var[NNESTS].i-1)*sizeof(double));
      x = opti.bf->x;

      grad = opti.bf->grad;

      amlet_init_bfgs(opti.bf, &dr);
 
      opti.bf->max_iter = niter;
      break;
    default:
      opti.b = NULL;
      c_error("unsupported method invoked for optimization process");
    }

    /* Alternatives sampling -- beta function */
    if (dr.var[ALTERNATIVES_SAMPLING].i) {
      printf("Alternatives sampling\n");
      amlet_alternatives_sampling(&dr, n_obs, n_availability,
				  n_choices, data.np,
				  data.nobs, data.nalt, data.nvar, n_rep,
				  dr.var[ALTERNATIVES_SAMPLING].i);
      printf("Alternatives sampling completed\n");
      exit(0);
    }

    /* Initialize data structure */
    data.obs = n_obs;
    data.rep = n_rep;
    data.availability = n_availability;
    data.choices = n_choices;

    /* Optimization process */
    elapsed_time = ctools_timer(mclog->f, &sec, &ms);

    if (!n_random_coeff) {
      if (m->verbose) printf("** Multinomial Logit **\n");
      nt_log_section(mclog, "Multinomial logit");

      switch (dr.var[METHOD].i) {
      case MCBTR:
      case MCBTRC:
      case MCBTRDA:
	ll = nlp_btrx(opti.b,
		      (C_GENERIC)wrapper_grad, NlpFG,
		      (C_GENERIC)wrapper_grad, NlpFval,
		      NULL, 0,
		      &data, 0.0,  mclog->f, m->I,
		      CblasRowMajor, CblasUpper, work);
	niter = opti.b->k;
	break;
      case MCRTR:
      case MCRTRDA:
	ll = nlp_rtrx(opti.b,
		      (C_GENERIC)wrapper_grad, NlpFG,
		      (C_GENERIC)wrapper_grad, NlpFval,
		      NULL, 0,
		      &data, 0.0,  mclog->f, m->I,
		      CblasRowMajor, CblasUpper, work);
	niter = opti.b->k;
	break;
      case MCBFGS:
      case MCBFGSDA:
	ll = nlp_bfgs(opti.bf, (C_GENERIC)wrapper, (C_GENERIC)wrapper_grad,
		      &data, mclog->f, m->I, work);
	niter = opti.bf->k;
	break;
      }
    }
    else {
      if (m->verbose) printf("** Mixed Logit **\n");
      nt_log_section(mclog, "Mixed logit");

      /* we construct an array that summarizes the basic types for random
	 variables */
      sampling_types = _mallocs(n_random_coeff, int);
      moca_set_sampling_types(btype, sampling_types, n_random_coeff);

      if (dr.var[NREP].i < MINSAMPLING) {
	c_warning("sampling size changed to the minimum allowed one (%d)",
		  MINSAMPLING);
	dr.var[NREP].i = MINSAMPLING;
      }

      /* memory allocation for (quasi-)random draws */
      data.sampling = amlet_sample_alloc(dr.var[QMCREP].i*data.np,
					 dr.var[NREP].i, n_random_coeff);
      data.size_sampling = data.max_size = dr.var[NREP].i;

      np = data.np;
      if (!dr.var[DRAWS].i) {
	/* Read external draws */
	f = c_fopen(dr.var[DRAWSNAME].s, "r");
	if (!f) {
	  c_warning("%s\n", amlet_info_message(AMLET_ERROR_DRAWS_FILE_OPENING));
	  dr.var[DRAWS].i = 1;
	}
	else {
	  if (ran_import_sample(f, np, dr.var[NREP].i*n_random_coeff,
				data.sampling)) {
	    c_warning("%s\n", amlet_info_message(AMLET_ERROR_DRAWS_FILE_READING));
	    dr.var[DRAWS].i = 1;
	  }
	  fclose(f);
	}
	if (dr.var[DRAWS].i > 0) {
	  printf("Overwrite %s with news draws? ([Y]/N)", dr.var[DRAWSNAME].s);
	  i = getchar();
	  if (i == 'N' || i == 'n') {
	    free(dr.var[DRAWSNAME].s);
	    dr.var[DRAWSNAME].s = 0;
	  }
	}
      }

      for (i = 0; i < 2; i++) ran_random_resetnextsubstream(ran);
      f = NULL;
      if (dr.var[DRAWS].i) {
	/* Generate new draws */
	if (dr.var[DRAWSNAME].s) {
	  f = fopen(dr.var[DRAWSNAME].s, "w");
	  if (!f) amlet_warning(mclog->f, 1,
				"unable to open %s to write draws\n",
				dr.var[DRAWSNAME].s);
	}
      }
      amlet_sample_generate(f, ran, &dr, dr.var[DRAWS].i,
			    data.sampling, sampling_types, n_random_coeff);
      if (f) fclose(f);

      /* Adjust configuration of optimizers */
      switch (dr.var[METHOD].i) {
      case MCBTR:
      case MCRTR:
      case MCBTRC:
      case MCRTRDA:
      case MCBTRDA:
	opti.b->d1 = wrapper_simulated_grad;
	if (m->verbose) amlet_log("%s\n", btr_get_step_name(opti.b));
	/*	  opti.b->verbose = HIGHLY_VERBOSE;  */
	break;
      case MCBFGS:
      case MCBFGSDA:
	opti.bf->d1 = wrapper_simulated_grad;
	break;
      }

      /* Initial sample size; this is set to a lower value if external
	 variable sampling size is used. */
      data.size_sampling =  dr.var[NREP].i;

      elapsed_time = ctools_timer(mclog->f, &sec, &ms);
      aes = amlet_external_strategy_new(dr.var[NREP].i,
					dr.var[NREP].i,
					dr.var[EPS].d);
      do {
	amlet_external_strategy_update(aes, &opti);

	/* Do some reinitialization tasks */
	switch (dr.var[METHOD].i) {
	case MCBTR:
	case MCRTR:
	case MCBTRC:
	case MCBTRDA:
	case MCRTRDA:
	  btr_set_delta(opti.b, AMLET_DEFAULT_TRUST_REGION_RADIUS);
	  break;
	default:
	  break;
	}
	
	tol = aes->tol;

	switch (dr.var[METHOD].i) {
	case MCBTR:
	  opti.b->tol = tol;
	  opti.b->output = stdout;
	  ll = nlp_btrx(opti.b,
			(C_GENERIC)wrapper_simulated_grad, NlpFG,
			(C_GENERIC)wrapper_simulated_grad, NlpFval,
			NULL, 0,
			&data, 0.0, mclog->f,
			m->I,
			CblasRowMajor, CblasUpper, work);
	  niter = opti.b->k;
	  break;
	case MCRTR:
	  opti.b->tol = tol;
	  opti.b->output = stdout;
	  ll = nlp_rtrx(opti.b,
			(C_GENERIC)wrapper_simulated_grad, NlpFG,
			(C_GENERIC)wrapper_simulated_grad, NlpFval,
			NULL, 0,
			&data, 0.0, mclog->f,
			m->I,
			CblasRowMajor, CblasUpper, work);
	  niter = opti.b->k;
	  break;
	case MCBTRDA:
	  opti.b->tol = tol;
	  /* Be careful: wrapper_accuracy is not really designed. */
	  ll = model_amlet_btr_da(wrapper_simulated_accuracy,
				  wrapper_simulated_accuracy_grad,
				  &data, 0.0, opti.b, m->I,
				  CblasRowMajor, CblasUpper,
				  work, &accuracy);

	  niter = opti.b->k;
	  break;
	case MCRTRDA:
	  opti.b->tol = tol;
	  /* Be careful: wrapper_accuracy is not really designed. */
	  ll = model_amlet_rtr_da(wrapper_simulated_accuracy,
				  wrapper_simulated_accuracy_grad,
				  &data, 0.0, opti.b, m->I,
				  CblasRowMajor, CblasUpper,
				  work, &accuracy);
	  niter = opti.b->k;
	  break;
	case MCBTRC: {
	  /* To clean and to merge!!! */
	  NLP_Constraint *c = NULL;
	  int *cidx;
	  int nidx = ds_list_length(data.np_list);
	  int t;
	  int l = 0;

	  opti.b->tol = tol;
	  opti.b->stop = btr_constrained_projection_stopping_test;
	  opti.b->step_function = nlp_spcgm_st;
	  opti.b->projection_wrapper = nlp_btr_proj_standard;

	  cidx = _mallocs(2*nidx, int);

	  j = 0;
	  for (i = 0; i < data.nvar; i++) {
	    t = btype[i];
	    if (amlet_is_nonparametric(t)) {
	      niter = moca_count_parameters(t, data.np_list->data);
	      cidx[2*l] = j;
	      cidx[2*l+1] = j+niter;
	      l++;
	      j += niter;
	    }
	    else
	      if (t != FIXED) j += moca_count_parameters(t, NULL);
	  }
	  niter = 0;

	  /* Set the monotonicity constraints */
	  c = nlp_constraint_new(NLP_MONOTONICITY, nidx, cidx, 1);
	  nlp_btr_constraints_register(opti.b, c);

	  opti.b->verbose = HIGHLY_VERBOSE;
	  /* opti.b->verbose = MODERATE_VERBOSE; */

	  ll = nlp_btrx(opti.b,
			(C_GENERIC)wrapper_simulated_grad, NlpFG,
			(C_GENERIC)wrapper_simulated_grad, NlpFval,
			NULL, 0,
			&data, 0.0, mclog->f,
			m->I,
			CblasRowMajor, CblasUpper, work);
	  niter = opti.b->k;
	  
	  free(cidx);

	  break;
	}
	case MCBFGS:
	case MCBFGSDA:
	  opti.bf->tol = tol;
	  ll = nlp_bfgs(opti.bf, (C_GENERIC)wrapper_simulated,
			(C_GENERIC)wrapper_simulated_grad, &data,
			mclog->f, m->I, work);
	  niter = opti.bf->k;
	  break;
	default:
	  c_error("unsupported method invoked for optimization process");
	  break;
	}
      } while (!amlet_external_strategy_stopping_test(aes));
      amlet_external_strategy_free(aes);
    }

    /* elapsed_time is set to the time consumed by the optimization process */
    elapsed_time = ctools_timer(mclog->f, &sec, &ms) - elapsed_time;
	  
    /* Post analysis */
    /* if statistics are required, compute the accuracy */
    if (accuracy < 0.0 && dr.var[STATISTICS].i && n_random_coeff) {
      if (m->verbose) printf("** Computing statistics **\n");
      sll_accuracy(x, &data, &accuracy, &bias, 0);
      if (m->verbose) printf("** End of computing statistics **\n");
    }
    /* bias of simulation */
    else if (accuracy >= 0.0)
      bias = simulation_bias(data.np, data.err_quantile, accuracy);

    /* Define the output file name */
    buff = set_output_name(m, buff, &dr, drivername);
    if (nsim > 1) {
      str = c_string_itoa(str, simulation);
      buff = c_string_append(buff, "_");
      buff = c_string_append(buff, str);
    }

    /* Compute ll0 and llc */
    special_likelihood_functions(n_random_coeff, &ll0, &llc,
				 &data, &dr, work);

    if (!dr.var[HESSIAN_APPROX].i)
      amlet_model_likelihood_hessian(m, x, grad, n_random_coeff, &data, work);
    else
      memcpy(*(m->H), *(m->I), p->n*p->n*sizeof(double));

    /* DANGER: ONLY THE UPPER PART IF FILLED AFTER THE OPTIMIZATION */

    /* Compute the information and the correlation matrices */
    nt_vect_zeros(2*p->n, work, 1);
    if (!n_random_coeff && data.np != data.nobs) {
      for (i = 0; i < data.np; i++) {
	memmove(data.sgrad[i], data.sgrad[n_rep[i]], p->n*sizeof(double));
	for (j = data.rep[i]+1; j < data.rep[i+1]; j++)
	  cblas_daxpy(p->n, 1.0, data.sgrad[j], 1, data.sgrad[i], 1);
      }
    } 

    /* compute the information matrix */
    bhhh_update(CblasRowMajor, CblasUpper, p->n, *(m->I),
		p->n, NULL, 1, NULL, 1, work);
    amlet_likelihood_Hinv(m);

    amlet_weak_covariance_matrix(m, p->n, 1.0/nsize);

    if (!n_random_coeff && data.np != data.nobs)
      amlet_robust_covariance_matrix(m, p->n, data.np*1.0/(nsize*nsize));
    else
      amlet_robust_covariance_matrix(m, p->n, 1.0/nsize);

    amlet_dump_information_matrices(m);

    dcl_combine_parameters(data.nvar, btype, x, data.np_list,
			   data.fixed, ff, work);

    if (solutions) op_solutions_add(solutions, ll, x);
 
    /* Print results */
    amlet_output(m, &dr, &ac, buff,
		 &data, work, ff,
		 dr.var[LEVEL].d, n_random_coeff,
		 ll, ll0, llc, cblas_dnrm2(npar, grad, 1),
		 seed, data.size_sampling,
		 accuracy, bias, elapsed_time, niter,
		 dr.var[METHOD].i);

    /* Results validation */
    if (dr.var[VALIDATION].i) {
      amlet_validate_results(&dr, npar, x, &data,
			     n_random_coeff, ran, sampling_types, 100);
    }
  }

  if (solutions) {
    f = fopen("sols.txt", "w");
    op_solutions_dump(solutions, f);
    fclose(f);
  }

 End:
  if (buff) free(buff);

  nt_log_section(mclog, "AMLET termination\n");
  nt_log_write(mclog, "%s\n", amlet_info_message(info));
  printf("\nAMLET termination: %s\n", amlet_info_message(info));

  if (m->verbose) printf("Free memory\n");

  /* Free memory */
  driver_clean(&dr);

  if (solutions) free(solutions);
  if (str) free(str);

  free(availability);
  free(n_availability);
  free(choices);
  free(n_choices);
  free(obs);
  free(n_obs);
  free(rep);
  free(n_rep);

  /* Free nests */
  amlet_nest_free(data.n);

  tmp = data.np_list;
  for (i = 0; i < data.nvar; i++) {  
    j = btype[i];
    if (amlet_is_nonparametric(j)) {
      amlet_nonparametric_free(tmp->data, j);
      tmp = tmp->next;
    }
  }
  ds_list_free(data.np_list, NULL);

  nt_log_free(mclog);
  if (ac.logname) free(ac.logname);

  ran_random_free(ran);

  free(data.sprob);
  free(data.som);

  free(work);

  free(btype);

  if (opti.b) {
    switch(dr.var[METHOD].i) {
    case MCBTR:
    case MCBTRC:
    case MCBTRDA:
      op_btr_free(opti.b);
      free(opti.b);
      break;
    case MCBFGS:
    case MCBFGSDA:
      nlp_bfgs_free(opti.bf);
    free(opti.bf);
    break;
    }
  }

  if (data.sampling) {
    amlet_sample_free(data.sampling, data.np);
    free(sampling_types);
  }

  if (data.sgrad) nt_matrix_free(data.sgrad);
  if (data.somgrad) nt_matrix_free(data.somgrad);

  if (data.newgrad) nt_matrix_free(data.newgrad);
  free(data.fixed);
  if (data.beta) free(data.beta);

  // BOOTSTRAP!  dcl_memory_deallocate(data.nobs, data.nalt, choices, obs, availability);
  free(ff);

  if (m->verbose) printf("Thanks for using %s.\n", PROGRAM);

  amlet_model_free(m);

  return info;
}
