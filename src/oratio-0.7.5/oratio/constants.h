/**
 * constants.h
 *
 * ORATIO Library
 * Copyright (C) 2002-2012 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#ifndef __ORATIO_NTOOLS_CONSTANTS_H__
#define __ORATIO_NTOOLS_CONSTANTS_H__

/**
 * Numerical constants
 */

/**
 * Various constants are usually defined in math.h, so we first load
 * them if they exist.
 */
#include <math.h>

#ifndef M_E
#define M_E         2.71828182845904523536028747135      /* e */
#endif

#ifndef M_EULER
#define M_EULER     0.5772156649015328606065
#endif

#ifndef M_LN2
#define M_LN2       0.69314718055994530941723212146      /* ln(2) */
#endif

#ifndef M_LN4
#define M_LN4       1.3862943611198906;                  /* ln(4) */
#endif

#ifndef M_LN10
#define M_LN10      2.30258509299404568401799145468      /* ln(10) */
#endif

#ifndef M_LNPI
#define M_LNPI      1.14472988584940017414342735135      /* ln(pi) */
#endif

#ifndef M_LN_SQPI_2
#define M_LN_SQPI_2 -0.1207822376352453; // Ln(Sqrt(Pi)/2)
#endif

#ifndef M_LOG2E
#define M_LOG2E     1.44269504088896340735992468100      /* log_2 (e) */
#endif

#ifndef M_LOG10E
#define M_LOG10E    0.43429448190325182765112891892      /* log_10 (e) */
#endif

#ifndef M_PI
#define M_PI        3.14159265358979323846264338328      /* pi */
#endif

#ifndef M_PI_2
#define M_PI_2      1.57079632679489661923132169164      /* pi/2 */
#endif

#ifndef M_PI_4
#define M_PI_4      0.78539816339744830961566084582      /* pi/4 */
#endif

#ifndef M_1_PI
#define M_1_PI      0.31830988618379067153776752675      /* 1/pi */
#endif

#ifndef M_2_PI
#define M_2_PI      0.63661977236758134307553505349      /* 2/pi */
#endif

#ifndef M_SQPI_2
#define M_SQPI_2    0.88622692545275801                  /* sqrt(pi)/2 */
#endif

#ifndef M_SQRT1_2
#define M_SQRT1_2   0.70710678118654752440084436210      /* sqrt(1/2) */
#endif

#ifndef M_SQRT2
#define M_SQRT2	    1.41421356237309504880               /* sqrt(2) */
#endif

#ifndef M_SQRT3
#define M_SQRT3     1.73205080756887729352744634151      /* sqrt(3) */
#endif

#ifndef M_SQRTPI
#define M_SQRTPI    1.77245385090551602729816748334      /* sqrt(pi) */
#endif

#ifndef M_SQRTPI_2
#define M_SQRTPI_2  0.88622692545275801;                 /* sqrt(pi)/2 */
#endif

#ifndef M_2_SQRTPI
#define M_2_SQRTPI  1.12837916709551257389615890312      /* 2/sqrt(pi) */
#endif

#ifndef M_TWO17                                         /* 2^17 */
#define M_TWO17     131072.0
#endif

#ifndef M_TWO53                                         /* 2^53 */
#define M_TWO53     9007199254740992.0
#endif

#endif  /* __ORATIO_NTOOLS_CONSTANTS_H__ */
