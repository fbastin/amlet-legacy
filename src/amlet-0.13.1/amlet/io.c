/**
 * io.c
 *
 * AMLET (Another Mixed Logit Estimation Tool)
 * Copyright (C) 2002-2010 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include "info.h"
#include "io.h"

#include <ctype.h>
#include <stdlib.h>
#include <string.h>

#include "amlet_analysis.h"
#include "credits.h"
#include "log.h"

#include <config.h>

#include <oratio/oratio.h>
#include <oratio/ctimer.h>
#include <oratio/stat.h>
#include <amlet/amlet.h>
#include <amlet/dcl.h>

#define MAX_COLUMNS 6

enum {
  RAW = 0,
  HTML,
  XML,
  LATEX,
  END_FORMAT
};

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
		      AMLET_Driver *dr, char *name)
{
  int len;

  if (!dr->output) {
    dr->output = strdup(name);
    dr->output = c_string_append(dr->output, ".output");
  }

  len = strlen(dr->output);
  buff = _reallocs(buff, len+1, char);
  strcpy(buff,  dr->output);

  return buff;
}

/**
 * Define the log file name.
 */
char *set_logname(char **buff, AMLET_Driver *dr, amlet_config *ac)
{
  int len;

  len = strlen(AMLET_DEFAULT_LOG);
  *buff = _mallocs(len+1, char);
  strcpy(*buff, AMLET_DEFAULT_LOG);

  return *buff;
}

int amlet_error(amlet_config *ac,
		FILE *f, int verbose, int error, char *msg, ...)
{
  va_list curseur;

  va_start(curseur, msg);

  amlet_warning_va(f, verbose, "ERROR", msg, curseur);

  va_end(curseur);

  return error;
}

int amlet_print_msg(FILE *f, const char *msg, va_list ap)
{
  char *buff;
  int n = 1024;

  /* Create the message */
  buff = _mallocs(n, char);
  while (vsnprintf(buff, n, msg, ap) == -1) {
    n <<= 1;
    buff = _reallocs(buff, n, char);
  }

  fprintf(f, "%s", buff);
  //  if (verbose && f != stdout && f != stderr) printf("%s", buff);

  n = strlen(buff);

  free(buff);

  return n;
}

int print_msg(FILE *f, const char *msg, va_list ap)
{
  char *buff;
  int n = 1024;

  /* Create the message */
  buff = _mallocs(n, char);
  while (vsnprintf(buff, n, msg, ap) == -1) {
    n <<= 1;
    buff = _reallocs(buff, n, char);
  }

  fprintf(f, "%s", buff);

  n = strlen(buff);

  free(buff);

  return n;
}

void line(FILE *output, int len)
{
  int i;

  for (i = 0; i < len; i++) fprintf(output, "-");
  fprintf(output, "\n");
}

void section(FILE *output, char *string)
{
  int i;

  fprintf(output, "\n%s\n", string);
  for (i = 0; i < strlen(string); i++)
    fprintf(output, "-");
  fprintf(output, "\n\n");
}

void subsection(FILE *output, char *string)
{
  int i;

  fprintf(output, "\n%s\n", string);
  for (i = 0; i < strlen(string); i++)
    fprintf(output, "-");
  fprintf(output, "\n\n");
}

void table_line_end(FILE *output, int format)
{
  switch(format) {
  case RAW:
    fprintf(output, "\n");
    break;
  case HTML:
    fprintf(output, "</tr>\n");
    break;
  default:
    break;
  }
}

void table_cell(FILE *output, int format, int len, char *s, ...)
{
  int i, l;
  va_list curseur;

  va_start(curseur, s);

  switch(format) {
  case RAW:
    l = print_msg(output, s, curseur);
    for (i = l; i < len; i++) fprintf(output, " ");
    break;
  case HTML:
    fprintf(output, "<td align=center>");
    if (strcmp(s, ""))
      print_msg(output, s, curseur);
    else
      print_msg(output, "&nbsp;", curseur);
    fprintf(output, "</td>\n");
    break;
  default:
    break;
  }

  va_end(curseur);
}

void table_header(FILE *output, int format, int len, char *s, ...)
{
  int i, l;
  va_list curseur;

  va_start(curseur, s);

  switch(format) {
  case RAW:
    l = print_msg(output, s, curseur);
    for (i = l; i < len; i++) fprintf(output, " ");
    break;
  case HTML:
    fprintf(output, "<th align=center><b>");
    if (strcmp(s, ""))
      print_msg(output, s, curseur);
    else
      print_msg(output, "&nbsp;", curseur);
    fprintf(output, "</b></th>\n");
    break;
  default:
    break;
  }

  va_end(curseur);
}

/**
 * Print on output the result of one iteration of the BTR algorithm
 */
void amlet_print_iteration_btr(FILE *output, BTR *btr, void *data)
{
  int i;

  if (output != stdout)
    amlet_print_iteration_btr(stdout, btr, data);
  fprintf(output, "\n");
  fprintf(output, "** Iteration %d **\n", btr->k);
  fprintf(output, "likelihood value: %17.15g\n", btr->value);
  fprintf(output, "theta: ");
  for (i = 0; i < btr->n; i++) fprintf(output, " %g", btr->x[i]);
  fprintf(output, "\n");
  fprintf(output, "norm of gradient: %g\n",
	  cblas_dnrm2(btr->n, btr->grad, 1));
  fprintf(output, "trust-region radius: %g\n", btr->delta);
  fprintf(output, "sampling size: %d\n",
	  ((dcl_data *)data)->size_sampling);
  /*  if (btr->verbose == HIGHLY_VERBOSE) {
    fprintf(output, "rho: %g\n", btr->rho);
    nt_vect_print(output, "s", btr->step, btr->n);
    nt_vect_print(output, "g", btr->grad, btr->n);
    } */
  fprintf(output, "\n");
}

/**
 * Print on output the result of one iteration of the BFGS algorithm
 */
void amlet_print_iteration_bfgs(FILE *output, BFGS *bf, void *data)
{
  int i;

  if (output != stdout)
    amlet_print_iteration_bfgs(stdout, bf, data);
  fprintf(output, "\n");
  fprintf(output, "** Iteration %d **\n", bf->k);
  fprintf(output, "theta: ");
  for (i = 0; i < bf->n; i++) fprintf(output, " %g", bf->x[i]);
  fprintf(output, "\n");
  fprintf(output, "norm of gradient: %g\n",
	  cblas_dnrm2(bf->n, bf->grad, 1));
  fprintf(output, "likelihood value: %g\n", bf->value);
  fprintf(output, "sampling size: %d\n",
	  ((dcl_data *)data)->size_sampling);
  fprintf(output, "\n");
}

/**
 * INPUT
 */

/**
 * Read vectors of draws stored in a file.
 *
 * @param f is the file where the draws have to be read.
 * @param n is the number of vectors to read.
 * @param s it the dimension of each vector to read.
 * @param S is the array where the draws will be saved.
 */
int moca_import_draws(FILE *f, int n, int s, double **S)
{
  int i, j;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      if (fscanf(f, "%lg", S[i]+j) != 1) return AMLET_ERROR_DRAWS_FILE_READING;

  return 0;
}

/**
 * Read previously generated (0,1) draws.
 *
 * @return 0 if the reading has proceeded correctly, an error code else.
 */
/*int amlet_input_draws(char *input, double **draws, int *type,
		      const int ndraws, const int ncoeff)
{
  int i, j;
  char *line, *buf, *end;

  double v;

  FILE *in;

  in = c_fopen(input, "r");
  if (!in) return AMLET_ERROR_DRAWS_FILE_OPENING;

  for (i = 0; i < ndraws; i++) {
    line = c_file_get_line(in);
    if (!line) return AMLET_ERROR_FILE_REALS;
    buf = line;
    for (j = 0; j < ncoeff; j++) {
      v = strtod(buf, &end);
      if (buf == end) {
	free(line);
	return AMLET_ERROR_FILE_REALS;
      }
      buf = end;
    }
    free(line);
  }

  fclose(in);

  return 0;
  } */

/**
 * OUTPUT
 */

void amlet_output_draws(char *output, double **draws,
			const int ndraws, const int ncoeff)
{
  FILE *out;

  unsigned int i, j;

  out = c_fopen(output, "w");
  if (!out) {
    c_warning("Unable to open %s to write draws; skip...\n");
    return;
  }

  for (i = 0; i < ndraws; i++) {
    for (j = 0; j < ncoeff; j++) fprintf(out, " %18.16g", draws[i][j]);
    fprintf(out, "\n");
  }

  fclose(out);
}

/**
 * Print the estimation results.
 *
 * @param m is the estimated model.
 *        output: name of the output file
 *        md: data structure of the model
 *        theta: estimated parameters
 * @param alpha: level of significance for statistical tests
 *        n_random_coeff: number of random coefficients
 *        ll: weighted likelihood value at the solution
 * @param ll0 is the likelihood value if all parameters are set to zero.
 *        llc: likelihood value if only constant parameters are considered
 *        gnorm: gradient norm at the solution
 *        seed: seed used for the random generator (only for mixed-logit)
 *        ndraws: number of random draws (only for mixed logit)
 *        accuracy: if positive, accuracy in the likelihood
 *                  evaluation, else not significant
 *        bias: if negative, bias of simulation
 *        elapsed_time: time consumed during the optimization process
 * @param niter is the number of iterations
 *        method: identifier of the used optimization algorithm
 *        format: format of ouput file
 */
void amlet_output_results(AMLET_Model *m, AMLET_Driver *dr,
			  FILE *f, dcl_data *md,
			  amlet_config *ac,
			  double *theta, int *fixed, double alpha,
			  int n_random_coeff,
			  double ll, double ll0, double llc,
			  double gnorm, long seed, int ndraws,
			  double accuracy, double bias,
			  int elapsed_time, int niter, int method,
			  int format)
{
  OP_Problem *problem = amlet_get_problem(m);

  int i, j;
  int len[10];
  double *t, *tmp;;
  double p;

  long sec, ms;

  /* Number of terms over which ll has been averaged: individuals for
     mixed logit models, observations for multinomial logit models.
     ll0 and llc are always averaged over observations. With panel data
     (np < nobs) the totals were previously all computed as mean*np. */
  int nll = n_random_coeff ? md->np : md->nobs;

  t = _mallocs(md->npar, double);
  tmp = _mallocs(md->npar, double);
  nt_vect_zeros(md->npar, tmp, 1);

  /* compute the t-statistics */
  p = amlet_t_statistics(md->npar, theta, tmp, m->RCov, alpha, t);

  free(tmp);

  /* Remember that we compute the opposite of the log-likelihood
     since we perform a minimization. */
  ll  = -ll;
  ll0 = -ll0;
  llc = -llc;

  amlet_credits(f, "", "");
  op_problem_print_name(f, problem);

  /* Problem characteristics */

  ctools_timer(f, &sec, &ms);

  section(f, "Problem specifications");

  len[1] = 10;
  len[2] = 15;
  len[3] = 15;
  len[0] = 0;
  for (i = 1; i <= 3; i++) len[0] += len[i];

  line(f, len[0]);
  table_header(f, format, len[1], "");
  table_header(f, format, len[2], "Individuals");
  table_header(f, format, len[3], "Observations");
  fprintf(f, "\n");
  line(f, len[0]);
  table_cell(f, format, len[1], "Number");
  table_cell(f, format, len[2], "%d", md->np);
  table_cell(f, format, len[3], "%d", md->nobs);
  table_line_end(f, format);
  /*  table_cell(f, format, len[1], "Rejected");
  table_cell(f, format, len[2], "%d", dr->var[NP].i-md->np);
  table_cell(f, format, len[3], "%d", dr->var[NOBS].i-md->nobs); */
  table_line_end(f, format);
  line(f, len[0]);

  fprintf(f, "\n");
  fprintf(f, "Number of random coefficients: %d\n", n_random_coeff);

  /* Information about execution */

  section(f, "Information about execution");

  fprintf(f, "Method used for the optimization: %s\n", methods[method].name);
  fprintf(f, "\n");

  len[1] = 10;
  len[2] = 23;
  len[3] = 12;
  len[0] = 0;
  for (i = 1; i <= 3; i++) len[0] += len[i];

  line(f, len[0]);
  table_header(f, format, len[1], "");
  table_header(f, format, len[2], "Optimization process");
  table_header(f, format, len[3], "Total");
  fprintf(f, "\n");
  line(f, len[0]);
  table_cell(f, format, len[1], "Time (s)");
  table_cell(f, format, len[2], "          %d", elapsed_time);
  table_cell(f, format, len[3], "%ld.%ld", sec, ms);
  table_line_end(f, format);
  line(f, len[0]);

  fprintf(f, "\n");
  fprintf(f, "Number of iterations: %d\n", niter);
  /*  fprintf(f, "Function evaluations: %d", nfctev);
      endofline(f, format); */

  if (n_random_coeff) {
    /* mixed-logit model */
    fprintf(f, "Draws generation method: %s\n",
	    moca_get_method_name(amlet_get_sampling_id(m)));
    fprintf(f, "Final number of random draws: %d\n", ndraws);
  }

  fprintf(f, "\n");

  /* Results */

  section(f, "Results and statistics");

  len[1] = 16;
  len[2] = 15;
  len[3] = 15;
  len[4] = 25;
  len[0] = 0;
  for (i = 1; i <= 4; i++) len[0] += len[i];

  line(f, len[0]);
  table_header(f, format, len[1], "Log-likelihood");
  table_header(f, format, len[2], "Total");
  table_header(f, format, len[3], "Mean");
  table_header(f, format, len[4], "/ number of observations");
  fprintf(f, "\n");
  line(f, len[0]);
  table_cell(f, format, len[1], "Optimal");
  table_cell(f, format, len[2], "%12.10g", ll*nll);
  table_cell(f, format, len[3], "%12.10g", ll);
  table_cell(f, format, len[4], "%12.10g", (ll*nll)/md->nobs);
  table_line_end(f, format);
  table_cell(f, format, len[1], "At zero");
  table_cell(f, format, len[2], "%12.10g", ll0*md->nobs);
  table_cell(f, format, len[3], "%12.10g", (ll0*md->nobs)/nll);
  table_cell(f, format, len[4], "%12.10g", ll0);
  table_line_end(f, format);
  table_cell(f, format, len[1], "At constants");
  table_cell(f, format, len[2], "%12.10g", llc*md->nobs);
  table_cell(f, format, len[3], "%12.10g", (llc*md->nobs)/nll);
  table_cell(f, format, len[4], "%12.10g", llc);
  table_line_end(f, format);
  line(f, len[0]);

  fprintf(f, "\n");
  fprintf(f, "Gradient norm: %12.10g\n", gnorm);

  if (md->npar != md->nvar) {
      /* if the accuracy has been estimated, we write it */
      if (accuracy >= 0.0) {
	  fprintf(f, "Accuracy of the likelihood function evaluation (mean value): %g\n", accuracy);
      }
      if (bias < 0.0) fprintf(f, "Bias of simulation: %g\n", bias);
  }

  fprintf(f, "\n");

  len[1] = 6;
  len[2] = 15;
  len[3] = 15;
  len[4] = 15;
  len[5] = 5;
  len[0] = 0;
  for (i = 1; i <= 5; i++) len[0] += len[i];

  fprintf(f, "\n");
  for (i = 0; i < md->npar; i++)
    fprintf(f, "%10.9g\n", theta[i]);
  fprintf(f, "\n");

  line(f, len[0]);
  table_header(f, format, len[1], "");
  table_header(f, format, len[2], "Estimation");
  table_header(f, format, len[3], "Standard deviation");
  table_header(f, format, len[4], "t-statistic");
  table_header(f, format, len[5], "");
  fprintf(f, "\n");
  line(f, len[0]);

  j = 0;
  for (i = 0; i < md->npar+md->nfixed; i++) {
    table_cell(f, format, len[1], "%4d.", i);
      table_cell(f, format, len[2], "%10.8g", theta[i]);
    if (!fixed[i]) {
      table_cell(f, format, len[3], "%10.8g", sqrt(m->RCov[j][j]));
      table_cell(f, format, len[4], "%10g", t[j]);
      /* Test if the parameter is statistically significant */
      if (t[j] >= p && !isnan(t[j]))
	table_cell(f, format, len[5], "*");
      else
	table_cell(f, format, len[5], "");
      j++;
    }
    table_line_end(f, format);
  }
  line(f, len[0]);

  subsection(f, "Statistics");

  ll *= nll;
  ll0 *= md->nobs;

  fprintf(f, "\n");
  fprintf(f, "rho-squared: %g\n", 1-ll/ll0);
  fprintf(f, "\n");
  fprintf(f, "adjusted rho-squared: %g\n",
	  1-(ll-md->npar)/ll0);
  fprintf(f, "\n");

  subsection(f, "Likelihood ratio test");

  fprintf(f, "Chi-square value for l(0): %g\n", st_likelihood_chi2(ll, ll0));
  fprintf(f, "Level of rejection: %g\n",
	  st_likelihood_chi2_level(md->npar, ll, ll0));

  /* Print the variance-covariance matrix, which corresponds to the
     upper part of the matrix I. */
  subsection(f, "Variance-covariance matrix");

  if (dr->var[ROBUST].i)
    fprintf(f, "Robust estimation");
  else
    fprintf(f, "Standard estimation");
  fprintf(f, "\n");

  for (i = 0; i < md->npar; i++) {
    for (j = 0; j < md->npar; j++) {
      if ((j % MAX_COLUMNS) == 0 && md->npar > MAX_COLUMNS) {
	fprintf(f, "Line %d (columns %d to %d):\n",
		i+1, j+1, nt_min(j+MAX_COLUMNS, md->npar));
      }
      if (i >= j)
	fprintf(f, " %12g", m->RCov[j][i]);
      else
	fprintf(f, " %12g", m->RCov[i][j]);
      if (((j+1) % MAX_COLUMNS) == 0 && j+1 < md->npar)
	fprintf(f, "\n");
    }
    fprintf(f, "\n");
  }

  /* Print the correlation matrix, which corresponds to the
     lower part of the matrix I. */
  subsection(f, "Correlation matrix");

  for (i = 0; i < md->npar; i++) {
    for (j = 0; j < md->npar; j++) {
      if ((j % MAX_COLUMNS) == 0 && md->npar > MAX_COLUMNS)
	fprintf(f, "Line %d (columns %d to %d):\n",
		i+1, j+1, nt_min(j+MAX_COLUMNS, md->npar));
      if (i == j)
	fprintf(f, " %12g", 1.0);
      else if (i > j)
	fprintf(f, " %12g", m->RCov[i][j]);
      else
	fprintf(f, " %12g", m->RCov[j][i]);
      if (((j+1) % MAX_COLUMNS) == 0 && j+1 < md->npar)
	fprintf(f, "\n");
    }
    fprintf(f, "\n");
  }
}

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
		 int elapsed_time, int niter, int method)
{
  FILE *f;

  if (output) {
    f = c_fopen(output, "w");
    if (!f) {
      c_warning("unable to open %s to write results; standard output will be used\n", output);
      f = stdout;
    }
  }
  else
    f = stdout;

  amlet_output_results(m, dr, f, md, ac, theta, fixed,
		       alpha, n_random_coeff,
		       ll, ll0, llc, gnorm,
		       seed, ndraws,
		       accuracy, bias, elapsed_time, niter,
		       method, 0);

  if (f != stdout) fclose(f);

  return 0;
}

/**
 * Read the repetitions file.
 *
 * @param repname is the name of the file containing the repetitions. The
 *   ith entry in the file represents the number of observations associated
 *   to individual i.
 * @param rep is a pointer to an integer array of size at least
 *   equal to np+1. On output rep[i+1] will the contain the number of
 *   observations associated to indivual i.
 * @param np is number of individuals.
 *
 * @return 0 is the lecture has encountered no difficulties, else an
 *   error identifier.
 */
int amlet_read_repetitions(char *repname, int *rep, int np)
{
  FILE *f;
  int i, r;

  i = 0;
  f = fopen(repname, "r");

  if (!f) return AMLET_ERROR_REPETITIONS_FILE_OPENING;

  while (fscanf(f, "%d", &r) > 0 && i < np) rep[++i] = r;
  fclose(f);

  return (i == np) ? 0 : AMLET_ERROR_REPETITIONS_FILE_READING;
}
