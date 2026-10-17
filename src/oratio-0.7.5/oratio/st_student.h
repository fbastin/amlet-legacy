/**
 * st_student.h
 *
 * ORATIO Library - statistical functions
 * Copyright (C) 2001-2012 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __ORATIO_STAT_STUDENT_H__
#define __ORATIO_STAT_STUDENT_H__

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Compute the Student probability distribution function with nu
   * degrees of freedom, at x.
   *
   * @param x is the point where we want to evaluate the PDF.
   * @param nu is the number of degrees of freedom; nu must be positive.
   *
   * @return the value of the pdf of a Student distribution with nu degrees of
   *   freedom, at x.
   */
  double st_student_pdf(double x, double nu);

  /**
   * Compute the Student cumulative distribution function with nu
   * degrees of freedom, at x.
   *
   * @param x is the point where we want to evaluate the CDF.
   * @param nu is the number of degrees of freedom; nu must be positive.
   *
   * @return the value of the cdf of a Student distribution with nu degrees of
   *   freedom, at x.
   */
  double st_student_cdf(double x, double nu);

  /**
   * Compute the Student inverse cumulative distribution function with nu
   * degrees of freedom, at u.
   *
   * @param u is the point where we want to evaluate the inverse CDF; u must be
   *   in [0,1].
   * @param nu is the number of degrees of freedom; nu must be positive.
   *
   * @return the value of the inverse cumulative distribution function of a
   *   Student with nu degrees of freedom, at u.
   */
  double st_student_icdf(double u, double nu);

#ifdef __cplusplus
}
#endif

#endif  /* __ORATIO_STAT_STUDENT_H__ */
