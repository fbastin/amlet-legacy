/**
 * io.h
 *
 * AMLET (Another Mixed Logit Estimation Tool)
 * Copyright (C) 2002-2009 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#ifndef __AMLET_IO_H__
#define __AMLET_IO_H__

#include <stdarg.h>
#include <stdio.h>

#include <ophelia/nlp.h>

#include "amlet.h"
#include "amlet_config.h"
#include "amlet_driver.h"

#ifdef __cplusplus
extern "C" {
#endif


  /**
   * OUTPUT
   */

  /**
   * Define the ouput file name.
   *
   * @param m is the model under consideration.
   * @param buff is a string pointer, that will replaced by a pointer to
   *   the newly formed string.
   */
  char *set_output_name(AMLET_Model *m, char *buff,
			AMLET_Driver *dr, char *name);

  /**
   * Define the log file name.
   */
  char *set_logname(char **buff, AMLET_Driver *dr, amlet_config *ac);

  int amlet_error(amlet_config *ac,
		  FILE *f, int verbose, int error, char *msg, ...);

  int amlet_print_msg(FILE *f, const char *msg, va_list ap);

  /**
   * Print on output the result of one iteration of the BTR algorithm.
   */
  void amlet_print_iteration_btr(FILE *output, BTR *btr, void *data);

  /**
   * Print on output the result of one iteration of the BFGS algorithm.
   */
  void amlet_print_iteration_bfgs(FILE *output, BFGS *bf, void *data);

  /**
   * Read previously generated (0,1) draws.
   *
   * @return 0 if the reading has proceeded correctly, an error code else.
   */
  //  int amlet_input_draws(char *input, double **draws, int *type,
  //			const int ndraws, const int ncoeff);

  void amlet_output_draws(char *output, double **draws,
			  const int nind, const int ncoeff);

  void amlet_output_header(FILE *output, int format);

  void amlet_print_credits(FILE *f, int format);

  /**
   * Results output.
   *
   * @param m is the estimated model.
   * @param ac is a pointet to the AMLET config structure; the content
   *   is unchanged on output.
   * @param output is the name of the output file.
   */
  int amlet_output(AMLET_Model *m, AMLET_Driver *dr,
		   amlet_config *ac, char *output,
		   dcl_data *md,
		   double *theta, int *fixed, double alpha,
		   int n_random_coeff,
		   double ll, double ll0, double llc,
		   double gnorm, long seed, int ndraws,
		   double accuracy, double bias,
		   int elapsed_time, int niter, int method);

#ifdef __cplusplus
}
#endif

#endif  /* __AMLET_IO_H__ */
