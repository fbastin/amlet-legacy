/**
 * nlp.h
 *
 * OPHELIA Library
 * Copyright (C) 1999-2008 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#ifndef __NLP_H__
#define __NLP_H__

/* COMMON routines */
#include "defaults.h"
#include "io.h"
#include "model.h"
#include "projection.h"
#include "stopping.h"

/* Simple methods */
#include "gss.h"
#include "sa.h"

/* Constraints */
#include "bounds.h"
#include "constraints.h"

/* BFGS algorithm */
#include "bfgs.h"
#include "linesearch.h"

/* BTR algorithm */
#include "btr.h"
#include "btr_struct.h"
#include "radius.h"
#include "pcg.h"
#include "pcg_monotonicity.h"
#include "step.h"

#include "rtr.h"

#include "outer.h"

/* Secant updates */
#include "secant.h"
#include "sr1.h"

/* Projections */
#include "isotonic.h"

/* Meta-heuristics */
#include "sa.h"

/* Miscellaneous */
#include "barrier.h"

#endif  /* __NLP_H__ */
