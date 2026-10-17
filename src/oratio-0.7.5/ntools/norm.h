/*****************************************************************************
 *
 * NTools Library
 * Copyright (C) 2002 Fabian Bastin <fbas@math.fundp.ac.be>
 *
 * Licensed under the Open Software License version 2.1
 *
 *****************************************************************************/

#ifndef __NTOOLS_NORM_H__
#define __NTOOLS_NORM_H__

/*****************************************************************************
 *
 * Norms computations
 *
 *****************************************************************************/

/* The 2-norm can be already calculated with a BLAS library, so we
   don't compute it here. We can however freely implement other
   norms. */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

  /**************************************************************************
   * 1-Norm
   **************************************************************************/
  double nt_vect_nrm1(int n, double *x, int inc);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __NTOOLS_NORM_H__ */
