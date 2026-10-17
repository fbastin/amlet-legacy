/**
 * io.c
 *
 * NLP (NonLinear Programming) Library
 * Copyright (C) 1999-2010 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

/**
 * Input/Ouput functions
 */

#include "io.h"
#include "stopping.h"
#include "problem.h"

/**
 * Print the problem name in the output file.
 *
 * @param output is the file where to print name.
 * @param p is the problem under consideration.
 */
void op_problem_print_name(FILE *output, OP_Problem *p)
{
  char *name = op_problem_get_name(p);

  if (name)
    fprintf(output, "* Problem %s *\n\n", name);
  else
    fprintf(output, "* Unamed problem *\n\n");
}


/**
 * Results
 */
int nlp_print_results(FILE *output, int format,
		      int methodid, void *method,
		      void *data)
{
  return 0;
}

/**
 * Iterations
 */

/**
 * General iteration output function
 *
 * #PRE:  n: problem dimension
 *        x: variables vector
 *        val: function value at x
 * @param k is the iteration index
 * @param gr is the gradient at the point x
 */
void general_print_iteration(FILE *output, int k, int n,
				    double *x, double *gr, double val)
{
  int i;

  fprintf(output, "** Iteration %d **\n", k);
  fprintf(output, "Objective: %12.10g\n", val);
  fprintf(output, "x:");
  for (i = 0; i <  n; i++) fprintf(output, " %g", x[i]);
  fprintf(output, "\n");  
  fprintf(output, "gradient:");
  for (i = 0; i <  n; i++) fprintf(output, " %g", gr[i]);
  fprintf(output, "\n");  
  fprintf(output, "gradient norm: %g\n", cblas_dnrm2(n, gr, 1));
}

/**
 * Outer iterations
 */

void nlp_outer_print_iteration(FILE *output, nlp_outer *no, void *data)
{
  fprintf(output, "\n");
  fprintf(output, "** Outer iteration %d **\n", no->k);
  fprintf(output, "Barrier parameter: %g\n", no->mu);
  fprintf(output, "Total of inner iterations: %d\n", no->t);
  fprintf(output, "\n");
}

/**
 * BFGS I/O functions
 */
void bfgs_print_iteration(FILE *output, BFGS *bf, void *data)
{
  general_print_iteration(output, bf->k, bf->n, bf->x, bf->grad, bf->value);
  fprintf(output, "\n");
}

void bfgs_print_results(FILE *output, BFGS *bf, void *data)
{
  int i;

  for (i = 0; i < 50; i++) fprintf(output, "-");
  fprintf(output, "\n");

  if (nlp_has_converged(bf->info))
    fprintf(output, "Convergence reached due to %s\n",
	    nlp_stopping_msg(bf->info));    
  else
    fprintf(output, "Convergence failed due to %s\n",
	    nlp_stopping_msg(bf->info));

  fprintf(output, "\nFinal value: %g\n", bf->value);
  fprintf(output, "Number of iterations: %d\n", bf->k);
  fprintf(output, "Final iterate:");
  for (i = 0; i < bf->n; i++) fprintf(output, " %g", bf->x[i]);
  fprintf(output, "\n");

  for (i = 0; i < 50; i++) fprintf(output, "-");
  fprintf(output, "\n");
}

/**
 * BTR I/O functions
 */

void btr_print_results(FILE *output, BTR *btr, void *data)
{
  int i;

  for (i = 0; i < 50; i++) fprintf(output, "-");
  fprintf(output, "\n");

  if (nlp_has_converged(btr->info))
    fprintf(output, "Convergence reached due to %s\n",
	    nlp_stopping_msg(btr->info));    
  else
    fprintf(output, "Convergence failed due to %s\n",
	    nlp_stopping_msg(btr->info));

  fprintf(output, "\nFinal value: %g\n", btr->value);
  fprintf(output, "Number of iterations: %d\n", btr->k);
  fprintf(output, "Final iterate:");
  for (i = 0; i < btr->n; i++)
    fprintf(output, " %g", btr->x[i]);
  fprintf(output, "\n");
  fprintf(output, "Trust-region radius: %10.8g\n", btr->delta);

  for (i = 0; i < 50; i++) fprintf(output, "-");
  fprintf(output, "\n");
}

void btr_print_iteration(FILE *output, BTR *btr, void *data)
{
  general_print_iteration(output, btr->k, btr->n, btr->x,
			  btr->grad, btr->value);
  fprintf(output, "Trust-Region Radius: %g\n", btr->delta);
  fprintf(output, "\n");
}

void btr_print_gnuplot(FILE *output, BTR *btr, void *data)
{
  int i;

  for (i = 0; i < btr->n; i++) fprintf(output, "%g ", btr->x[i]);
  fprintf(output, "%g\n", btr->value);
}

void btr_print_iteration_latex(FILE *output, BTR *btr, void *data)
{
  int i;

  fprintf(output, "%d &", btr->k);
  for (i = 0; i < btr->n; i++) fprintf(output, " %g &", btr->x[i]);
  fprintf(output, " %g & %g \\\\\n", btr->delta, btr->value);
  fprintf(output, "\\hline \n");
}
