/**
 * st_fisher.h
 *
 * ORATIO Library - statistical module
 * Copyright (C) 2001-2012 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#ifndef __ORATIO_STAT_FISHER_H__
#define __ORATIO_STAT_FISHER_H__

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Fisher-Snedecor distribution.
   */

  /**
   * Evaluation of the Fisher-Snedecor function with parameters a and b,
   * at point x.
   */
  double st_fisher_pdf(double x, double a, double b);

  /**
   * Moyenne de la fonction de Fisher-Snedecor de paramètres a, b
   *
   * #PRE: b > 2
   */
  double st_fisher_mean(double x, double a, double b);

  /**
   * Variance de la fonction de Fisher-Snedecor de paramètres a, b
   *
   * #PRE: b > 4
   */
  double st_fisher_variance(double x, double a, double b);

#ifdef __cplusplus
}
#endif

#endif /* __ORATIO_STAT_FISHER_H__ */
