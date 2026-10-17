/**
 * validation.c
 *
 * AMLET (Another Mixed Logit Estimation Tool)
 * Copyright (C) 2003-2010 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include "validation.h"

#include <stdio.h>

#include <oratio/oratio.h>

#include "amlet.h"

/**
 * Validate estimation results.
 *
 * @param x is the estimated parameters vectors.
 * @param nrand is the number of random coefficients in the model.
 * @param ran is an initialized random numbers generator.
 * @param nsim is the number of likelihood simulations to perform.
 */
void amlet_validate_results(AMLET_Driver *dr, int n, double *x, dcl_data *md,
			    int nrand, Random *ran, int *sampling_types,
			    int nsim)
{
  double *lvalues, *lbias;
  double err, variance;

  int bias = md->bias;

  /* double epsilon */
  double meanl, meanb, stdv;

  FILE *f;

  int i;

  nsim = 200;
	
  printf("** Validation **\n");
  printf("Please be patient: this may take a while.\n");

  if (dr->output) {
    f = fopen(dr->output, "a");
    if (!f) {
      c_warning("Unable to open %s to write validation results;\nscreen output will be used\n");
      f = stdout;
    }
  }
  else
    f = stdout;

  fprintf(f, "\n");
  fprintf(f, "Validation of results\n");
  fprintf(f, "---------------------\n");
  fprintf(f, "\n");

  /* memory allocation */
  lvalues = _mallocs(nsim, double);
  lbias = _mallocs(nsim, double);

  double v;
	
  //  md->bias = 1;
  /* compute likelihood and gradient for different samplings */
  for (i = 0; i < nsim; i++) {
    amlet_sample_generate(NULL, ran, dr, 1, md->sampling,
			  sampling_types, nrand);
    /* computation of likelihood and gradient */
    md->bias = 0;
    md->antithetics = 0;
    lvalues[i] = -simulated_likelihood_accuracy
      (x, md->obs, md->choices, md->sampling, md->btype,
       md->size_sampling, md->size_sampling, md, &err, &variance, 0);
    lbias[i] = simulation_bias(md->np, md->err_quantile, err); 
    md->bias = 1;
    v = -simulated_likelihood_accuracy
      (x, md->obs, md->choices, md->sampling, md->btype,
       md->size_sampling, md->size_sampling, md, &err, &variance, 0);
    //    md->antithetics = 1;
	  
    fprintf(f, "%3d %18.16g %18.16g %18.16g\n", i+1, v, lvalues[i], lbias[i]);
  }
	md->bias = bias;

  st_mean_variance(lvalues, 1, 0, nsim, &meanl, &stdv);
  stdv = sqrt(stdv);

  fprintf(f, "\n");
  fprintf(f, "Mean likelihood value: %17.15g\n", meanl);
  fprintf(f, "Standard deviation: %17.15g\n", stdv);
  fprintf(f, "Confidence radius: %17.15g\n", md->err_quantile*stdv);

  st_mean_variance(lbias, 1, 0, nsim, &meanb, &stdv);
  fprintf(f, "Mean lbias value: %17.15g\n", meanb);
  fprintf(f, "Standard deviation: %17.15g\n", stdv);

  fprintf(f, "Covariance likelihood-bias: %17.15g\n",
	  st_covariance(nsim, lvalues, 1, meanl, lbias, 1, meanb, 0, nsim));

  if (f != stdout) fclose(f);

  /* free memory */
  free(lvalues);
  free(lbias);

  printf("End of validation process.\n");
}
