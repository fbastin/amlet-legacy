/**
 * amlet_btr.h
 *
 * AMLET (Another Mixed Logit Estimation Tool)
 * Copyright (C) 2002-2004 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#ifndef __AMLET_BTR_H__
#define __AMLET_BTR_H__

#include <amlet/data.h>
#include <amlet/dcl.h>
#include <ophelia/nlp.h>

#define MINSAMPLING 32  /* minimum sampling size; 2^5 */

#ifdef __cplusplus
extern "C" {
#endif

  /* Constant names of variable sampling sizes internal strategies */
  enum {
    MCBTRDA_MODEL = 0,
    MCBTRDA_FCT
  };

  void dcl_data_clone_btrda(dcl_data *dest, dcl_data *src);

  int amlet_btr_radius_update(BTR *b, double rho,
			      double *x_cand, double new_value);

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
			    double *work, double *epsilon);


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
			    double *work, double *epsilon);
#ifdef __cplusplus
}
#endif

#endif  /* __AMLET_BTR_H__ */
