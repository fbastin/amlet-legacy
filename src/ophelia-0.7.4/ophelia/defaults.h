/**
 * defaults.h
 *
 * NLP (NonLinear Programming) Library
 * Copyright (C) 1999-2005 Fabian Bastin <fbas@math.fundp.ac.be>
 *
 * Licensed under the MIT License
 */

#ifndef __NLP2_DEFAULTS_H__
#define __NLP2_DEFAULTS_H__

#define NLP_MAX_ITER      10000
#define NLP_TOL           1e-6
#define NLP_MINIMUM_TOL   1e-12   /* minium allowed tol */

/**
 * Default values for BTR parameters
 */

#define BTR_ETA_1         0.2
#define BTR_ETA_2         0.75
#define BTR_GAMMA_1       0.5
#define BTR_GAMMA_2       0.5

#define BTR_INITIAL_DELTA 1.0

#endif  /* __NLP2_DEFAULTS_H__ */
