/**
 * likelihood.h
 *
 * AMLET Library
 * Copyright (C) 2002-2012 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#ifndef __AMLET_LIKELIHOOD_H__
#define __AMLET_LIKELIHOOD_H__

#include <ophelia/data.h>

#ifdef __cplusplus
extern "C" {
#endif

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
		    int shift);

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
		    double *accuracy, double *bias, int start);

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
			      int ssamp, int max_size, dcl_data *md);

  /**
   * @param theta is the vector of model parameters.
   * @param md is the data structure about the problem.
   *
   * @return 0.0.
   */
  double likelihood_gradient(double *theta, int n, void *data, double *grad);

  /**
   * Likehilood functions with accuracy estimation
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
	       double *grad, int start);

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
				       double *variance,
				       int start);

  /**
   * Compute the bias of simulation
   *
   * @param popsize is the population size
   * @param accuracy is the estimated precision.
   * @param alpha is the significance level used for the accuracy estimation.
   *
   * @return the bias of simulation
   */
  double simulation_bias(int popsize, double alpha, double accuracy);

  /**
   * Wrappers
   */

  /**
   * @param theta is the vector of model parameters.
   * @param n is the size of theta, i.e. the number of parameters.
   * @param data is a pointer to the data structure, of type amlet_data,
   *   describing the problem under investigation.
   */
  double wrapper(double *theta, int n, void *data);

  /**
   * @param theta is the vector of model parameters.
   * @param n is the size of theta, i.e. the number of parameters.
   * @param data is a pointer to the data structure, of type amlet_data,
   *   describing the problem under investigation.
   */
  double wrapper_grad(double *theta, int n, void *data, double *grad);

  /**
   * Compute the mixed-logit function.
   *
   * @param theta is the vector of model parameters.
   * @param n is the size of theta, i.e. the number of parameters.
   * @param data is a pointer to the data structure, of type amlet_data,
   *   describing the problem under investigation.
   */
  double wrapper_simulated(double *theta, int n, void *data);

  /**
   * Compute the mixed-logit function and its gradient.
   *
   * @param theta is the vector of model parameters.
   * @param n is the size of theta, i.e. the number of parameters.
   * @param data is a pointer to the data structure, of type amlet_data,
   *   describing the problem under investigation.
   */
  double wrapper_simulated_grad(double *theta, int n, void *data,
				double *grad);

  /**
   * @param theta is the vector of model parameters.
   * @param n is the size of theta, i.e. the number of parameters.
   * @param data is a pointer to the data structure, of type amlet_data,
   *   describing the problem under investigation.
   */
  double wrapper_simulated_accuracy(double *theta, int n, void *data,
				    double *epsilon, double *variance,
				    int start);

  /**
   * @param theta is the vector of model parameters.
   * @param n is the size of theta, i.e. the number of parameters.
   * @param data is a pointer to the data structure, of type amlet_data,
   *   describing the problem under investigation.
   */
  double wrapper_simulated_accuracy_grad(double *theta, int n, void *data,
					 double *epsilon, double *variance,
					 double *grad, int start);

#ifdef __cplusplus
}
#endif

#endif  /* __AMLET_LIKELIHOOD_H__ */
