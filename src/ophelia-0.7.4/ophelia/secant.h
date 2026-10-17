/**
 * secant.h
 *
 * NLP (NonLinear Programming) Library
 * Copyright (C) 2005-2006 Fabian Bastin <bastin@cerfacs.fr>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __NLP2_SECANT_H__
#define __NLP2_SECANT_H__

#include <oratio/cblas.h>

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Secant update prototype function.
   *
   * @param Order designs the orientation of the matrix Bk. The
   *   possible values are CblasRowMajor and CblasColMajor.
   * @param Uplo designs the triangular part that has to be updated. The
   *   possible value CblasUpper and CblasLower.
   * @param n is the number of variables.
   *
   * @return 1 if the approximation Bk has been updated, 0 else.
   */
  typedef int (*nlp_secant)(const enum CBLAS_ORDER Order,
			    const enum CBLAS_UPLO Uplo,
			    int n, int m, double *hessian,
			    double *values, double *valuesplus,
			    double *grad, double *gradplus,
			    double *step,
			    int *info, double *work, void *data);

  enum {
    NLP_SECANT_NONE = 0,
    NLP_SECANT_BFGS,
    NLP_SECANT_SR1,
    NLP_GAUSS_NEWTON,
    NLP_SECANT_LS
  };

  typedef struct _NlpSecant NlpSecant;

  struct _NlpSecant {
    double **H;         /* approximation matrix */
    int n;              /* main dimension of H */
    int id;             /* approximation method identifier */
    nlp_secant update;  /* approximation update function */
    nlp_secant init;    /* approximation initialization function */
    int initialized;    /* 1 if the approximation is initialized, 0 else */
    NlpSecant *ns;      /* this may be useful for composite secante updates */
    void *data;
  };

  void *nlp_set_hessian_approximation(int id, int info);

  /**
   * @param n is the number of variables.
   */
  NlpSecant *nlp_secant_new_default(int n, int id);

  /**
   * @param n is the number of variables.
   */
  NlpSecant *nlp_secant_new(int n, int id);

  void nlp_secant_free(NlpSecant *ns);

  int nlp_secant_set_fct(NlpSecant *ns);

  /**
   * @param n is the number of variables.
   */
  int nlp_secant_init_default(const enum CBLAS_ORDER Order,
			      const enum CBLAS_UPLO Uplo,
			      int n, int m, double *H,
			      double *values, double *valueplus,
			      double *grad, double *gradplus,
			      double *step, int *info, double *work,
			      void *data);

  /**
   * @param n is the number of variables.
   */
  void nlp_secant_init(NlpSecant *ns, const enum CBLAS_ORDER Order,
		       const enum CBLAS_UPLO Uplo,
		       int n, int m,
		       double *grad, double *gradplus,
		       double *step, int *info, double *work, void *data);

  /**
   * Define the Hessian matrix associated to the particular secant method.
   *
   * @param n is the number of variables.
   *
   * @return a pointer to the Hessian matrix corresponding to the
   *   secant method defined in ns.
   */
  double **nlp_secant_get_hessian(NlpSecant *ns, const enum CBLAS_ORDER Order,
				  const enum CBLAS_UPLO Uplo, int n, int m);

#ifdef __cplusplus
}
#endif

#endif  /* __NLP2_SECANT_H__ */
