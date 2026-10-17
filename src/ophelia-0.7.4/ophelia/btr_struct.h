/**
 * btr_struct.h
 *
 * OPHELIA Library
 * Copyright (C) 1999-2008  Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#ifndef __OPHELIA_BTR_STRUCT_H__
#define __OPHELIA_BTR_STRUCT_H__

#include <stdio.h>
#include <sys/types.h>

#include <oratio/oratio.h>

#include "common_struct.h"
#include "constraints.h"
#include "projection.h"
#include "secant.h"

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Structure definition
   */

  typedef struct _BTR BTR;

  /**
   * Step computation
   */
  typedef void (*BTR_step)(BTR *b,
			   const enum CBLAS_ORDER Order,
			   const enum CBLAS_UPLO Uplo,
			   double **hessian, double *work,
			   void *data);

  /**
   * Projection functions.
   */
  typedef int (*BTR_proj)(BTR *btr, double *x);
  typedef double (*BTR_proj_norm)(BTR *btr, double *x);

  /**
   * Basic trust-region structure.
   *
   * The memory gestion of data, proj_data, prob_info is at the user
   * discretion.
   */
  struct _BTR {
    /* BTR algorithm constants */
    double eta1, eta2, gamma1, gamma2;

    int n;             /* problem dimension */
    int m;             /* second dimension, set to 0 if unused */
    double *x;         /* point */
    int *x_type;       /* composant types */

    double *scaling;

    double value;      /* function value */
    double *grad;      /* function gradient at x */

    double delta;      /* trust-region radius */
    void (*update_radius)(BTR *b, double rho);

    int retro;
    double rho;

    /* step towards next candidate */
    int step_type;     /* step computation identifier */
    double *step;
    BTR_step step_function;

    /* Iterations */
    int k;         /* iteration counter */
    int max_iter;  /* maximal number of iterations */

    /* Constraints */
    DS_List *constraints;  /* Constraints sets. */
    int nconstraints;      /* Number of constraints sets. */

    struct _NlpSecant **secant_updates;
    int n_secant_updates;

    /* derivation functions */
    void (*d2approx)(void); /* generic function for Hessian approximation - TEMPORARY !!! */

    /* derivation functions */
    deriv1 d1;     /* first derivate - TEMPORARY!!!!*/

    double (*wrapper_fct)(C_GENERIC f, int f_type,
			  double *x, int n, int m, void *data,
			  double *values, double *grad, double **hessian,
			  void *prob_info, int type_prob, int *info);
    void (*wrapper_grad)(C_GENERIC f, int f_type, C_GENERIC g, int g_type,
			 double *x, int n, void *data, double val,
			 double *grad, double *scaling,
			 void *prob_info, int prob_type, int *info,
			 double *h, double *work);
    void (*wrapper_hess)(const enum CBLAS_ORDER Order,
			 const enum CBLAS_UPLO Uplo,
			 C_GENERIC f, int f_type,
			 C_GENERIC g, int g_type,
			 C_GENERIC H, int h_type,
			 double *x , int n, void *data,
			 double *grad, double **hessian, double *scaling,
			 void *prob_info, int prob_type, int *info,
			 double *work);

    /* Projection functions, set to NULL if no projection has to
       be performed. */
    BTR_proj      projection_wrapper;
    BTR_proj_norm projection_norm_wrapper;
    /* Projection data; NULL if no projection has to be performed. */
    void*         projection_data;

    /* feasibility test */
    int (*feasible)(double *x, int n, void *data);

    /* Output functions */
    FILE *output;
    void (*printer)(FILE *output, BTR *btr, void *data);

    /* Stopping test */
    int (*stop)(BTR *btr, double rho);
    double tol;     /* tolerance */

    /* Trust-region radius initialization */
    double (*radius_init)(BTR *btr, void *data);

    /* info on the behaviour of the BTR algorithm */
    int info;

    void *data;
    void *prob_info;

    /* Distance to the boundary of the feasible set. */
    double (*dist)(double *x, int n, void *data);

    int verbose;

    double *work;   /* working vector */
  };

  /**
   * Initialization functions
   */

  /**
   * Initialize a BTR structure.
   *
   * @param btr is a memory allocated BTR structure.
   * @param n is the problem dimension.
   * @param m is the secondary, problem dependent, dimension, usually set to 0.
   */
  void op_btr_init(BTR *btr, int n, int m);

  /**
   * Reset a btr structure.
   *
   * @param btr is an initialized BTR structure.
   */
  void op_btr_reset(BTR *btr);

  /**
   * Create a new BTR structure.
   */
  BTR *op_btr_new();

  /**
   * Free memory associated to a BTR object.
   *
   * @param btr is an initialized BTR structure.
   */
  void op_btr_free(BTR *btr);

  /**
   * Setters
   */

  /**
   * Setting functions
   */

  void btr_set_constants(BTR *btr, double eta1, double eta2,
			 double gamma1, double gamma2);

  /**
   * Set the trust-region radius.
   *
   * @param btr is an initialized BTR structure.
   * @param delta is the desired value for the trust-region radius.
   */
  void btr_set_delta(BTR *btr, double delta);

  /**
   * Set the tolerance used to test convergence.
   */
  void btr_set_tol(BTR *btr, double tol);

  /**
   * Set the projection information.
   *
   * @param btr is an initialized BTR structure.
   */
  void btr_set_projection_data(BTR *btr, void *data);

  /**
   * Set the maximum iterations number.
   */
  void btr_set_max_iter(BTR *btr, int max_iter);

  void btr_set_hessian_approx(BTR* btr, C_GENERIC fct);

  void nlp_btr_wrapper_default(BTR *btr);

  /**
   * Secant updates management.
   */

  /**
   * @return the secant hessian approximation.
   */
  double **nlp_btr_secant_register(BTR *btr, int secant_id);

  double **btr_set_hessian_approximation(BTR *b, int id, int info);

  /**
   * Register a constraints block.
   *
   * @param btr is an initialized BTR structure.
   *
   * @return 1 if the registration has been successful 0 else.
   */
  int nlp_btr_constraints_register(BTR *btr, NLP_Constraint *c);

  /**
   * @param btr is an initialized BTR structure.
   */
  void nlp_btr_constraints_deregister_all(BTR *btr);

  /**
   * @param b is an initialized BTR structure.
   * @param y is the point to project.
   */
  int nlp_btr_proj_standard(BTR *b, double *y);

#ifdef __cplusplus
}
#endif

#endif  /* __OPHELIA_BTR_STRUCT_H__ */
