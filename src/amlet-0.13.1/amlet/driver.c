/**
 * driver.c
 *
 * AMLET (Another Mixed Logit Estimation Tool)
 * Copyright (C) 2002-2015 Fabian Bastin <fabian.bastin@gmail.com>
 *
 * Licensed under the MIT License
 */

/**
 * BUG TO FIX: censor only reads integers, not reals.
 */

#include "amlet.h"
#include "amlet_config.h"

#include "amlet_driver.h"
#include "driver_parser.h"

#include "amlet_data.h"
#include "nest.h"

#include "model.h"

#include <oratio/ctools.h>

#define DRIVER_DEFAULT_NREP 1000

/**
 * External parsing functions, defined in the Yacc source code driver_parser.y
 */

extern void amlet_parser_init(AMLET_Driver *d);
extern void amlet_parser_clean();
extern int yyparse ();

void driver_init_defaults(AMLET_Driver *dr)
{
  /**
   * Give default values for dr
   */

  int i;

  dr->var[NP].i = -1;
  dr->var[NALT].i = -1;
  dr->var[NOBS].i = 5000;
  dr->var[NVAR].i = 3;
  dr->var[IDDEP].i = 0;

  dr->var[CENSOR].i = 0;
  dr->var[IDCENSOR].i = 0;

  dr->var[RESCALE].i = 0;

  dr->var[BIAS].i = 1;

  dr->var[NITER].i = 1000;
  dr->var[EPS].d = AMLET_DEFAULT_TOL;
  dr->var[SEED1].i = 12345;
  dr->var[NREP].i = 0;
  dr->var[VARREP].i = 0;
  dr->var[LEVEL].d = 0.05;
  dr->var[LEVEL_ACCURACY].d = 0.10;

  dr->var[FORMAT].i = TRAIN96;

  dr->var[METHOD].i = MCBTRDA;
  dr->var[HESSIAN_UPDATE].i = HESSIAN_BFGS;
  dr->var[HESSIAN_APPROX].i = 0;
  dr->var[ROBUST].i = 1;

  dr->var[RANGEN].i = RAN_MRG32K3A;
  dr->var[LEGACY_RNG].i = 0;
  dr->var[SAMPLING].i = MOCA_MONTE_CARLO;
  dr->var[RANDOMIZATION].i = 0;
  dr->var[QMCREP].i = 1;

  dr->var[VALIDATION].i = 0;
  dr->var[STATISTICS].i = 0;

  dr->var[DRAWS].i = 1;

  dr->var[SYNTHETIC].i = 0;

  dr->var[ALTERNATIVES_SAMPLING].i = 0;

  for (i = 0; i < END_ARRAY_WORDS; i++) {
    dr->array[i].array = 0;
    dr->array[i].n = 0;
  }

  dr->generation = 0;
  dr->obsname = 0;
  dr->output = 0;

  dr->alt_types = 0;
  dr->alt_values = 0;

  /* nests */
  dr->var[NNESTS].i = 1;
  dr->nests = NULL;

  /* strings */
  for (i = NAME; i < END_RESERVED_WORDS; i++) dr->var[i].s = 0;
}

void driver_free_array(AMLET_Driver *dr, int id)
{
  free(dr->array[id].array);
  dr->array[id].array = NULL;
}

void driver_clean(AMLET_Driver *dr)
{
  int i;
  int nalt;

  nalt = dr->var[NALT].i;

  for (i = 0; i < END_ARRAY_WORDS; i++) driver_free_array(dr, i);
  if (dr->nests) {
    free(dr->nests);
  }

  if (dr->obsname) {
    free(dr->obsname);
    dr->obsname = 0;
  }
  if (dr->output) {
    free(dr->output);
    dr->output = 0;
  }

  if (dr->alt_types) {
    for (i = 0; i < nalt; i++) free(dr->alt_types[i]);
    for (i = 0; i < nalt; i++) free(dr->alt_values[i]);
    free(dr->alt_types);
    free(dr->alt_values);
  }
}

/**
 * CHECKING
 */

/**
 * Check the parameters defined in the driver
 *
 * @param dr is a pointer to the driver structure.
 * @param data is a pointer to the structure containing important data about
 *   the problem under consideration.
 */
void driver_check_parameters(AMLET_Driver *dr, amlet_config *ac, dcl_data *data)
{
  int i;

  if (!data->n_random_coeff) {
#ifdef DEBUG
    printf("** Multinomial model **\n");
#endif
    if (dr->var[VALIDATION].i) {
      c_warning("validation is designed for mixed-logit models only\n");
      dr->var[VALIDATION].i = 0;
    }
  }
  else {
#ifdef DEBUG
    printf("** Mixed logit model **\n");
#endif

    /* Check data format */
    if (dr->var[FORMAT].i <= 0 || dr->var[FORMAT].i >= END_DATA_FORMAT) {
      c_error("invalid data format");
    }
    
    /* Check parameters value */
    if (dr->var[NALT].i <= 1) {
      c_error("there must be at least two alternatives in the model");
    }

    /* Check number of draws per individual */
    if (dr->var[NREP].i < 10) {
      if (loose) {
	amlet_warning(stderr, 1,
		      "at least 10 draws per individual must be used; nrep set to %d", DRIVER_DEFAULT_NREP);
	dr->var[NREP].i = DRIVER_DEFAULT_NREP;
      }
      else {
	c_error("at least 10 draws per individual must be used");
      }
    }
  }

  if (dr->generation) {
    /* Test alt_types and alt_values */
    /* This check is still limited */
    if (!dr->alt_types || !dr->alt_values) {
      c_error("incomplete alternatives definitions");
    }
    else {
      for (i = 0; i <  dr->var[NALT].i; i++) {
	if (!dr->alt_types[i] || !dr->alt_values[i]) {
	  c_error("missing definition for alternative %d\n", i);
	}
      }
    }
  }
}

int driver_check_array_pos(AMLET_Driver *dr, int id, int pos)
{
  if (id >= ALTTYPES && id <= ALTVALUES) {
    if (!dr->generation) return 0;  /* should be more explicit */
    if (pos < 0 || pos >= dr->var[NALT].i) return 0;
  }
  else if (id == RESCALE) {
    if (pos < 0 || pos > 1) return 0;
  }

  return 1;
}

/**
 * Do some checks after the parsing
 */
void driver_check(AMLET_Driver *dr) {
  if (!dr->var[CENSOR].i) {
    if (dr->var[IDCENSOR].i) {
      c_warning("IDCENSOR set while CENSOR equal to 0; IDCENSOR reset to 0\n");
      dr->var[IDCENSOR].i = 0;
    }
  }
  if (dr->var[NITER].i < 0) {
    c_warning("number of iterations set to a negative value; reset to default value of %d\n", AMLET_DEFAULT_NITER);
    dr->var[NITER].i = AMLET_DEFAULT_NITER;
  }
  if (dr->var[EPS].d <= 0.0) {
    c_warning("negative or null tolerance; reset to default value of %g\n",
	      AMLET_DEFAULT_TOL);
    dr->var[EPS].d = AMLET_DEFAULT_TOL;
  }

  if (dr->generation) {
    if (!dr->alt_types || !dr->alt_values) {
      c_error("generation mode, but undefined alttypes or altvalues\n");
    }
  }

  /* If dynamic accuracy is used, we should not try to do external
     varying sampling size strategies */
  switch (dr->var[METHOD].i) {
  case MCBTRDA:
    if (dr->var[VARREP].i) {
      c_warning("external sampling size strategies are not compatible with dynamic accuracy procedures; varrep is reset to false.\n");
      dr->var[VARREP].i = 0;
    }
    break;
  }
}

int driver_set_string(AMLET_Driver *dr, int id, char *str)
{
  /* Be careful: could provocate a Segmentation Fault */
  if (dr->var[id].s)
    free(dr->var[id].s);
  else
    dr->var[id].s = str;

  return 1;
}

void driver_set_generation_process(AMLET_Driver *dr, char *obsname, int nalt)
{
  int i;

  if (dr->generation) {
    c_warning("generation process redefined\n");
    for (i = 0; i < nalt; i++) free(dr->alt_types[i]);
    for (i = 0; i < nalt; i++) free(dr->alt_values[i]);
    free(dr->alt_types);
    free(dr->alt_values);
  }

  dr->obsname = obsname;
  dr->alt_types = _mallocs(nalt, int *);
  dr->alt_values = _mallocs(nalt, double *);
}

/**
 * @param x is a double vector of size at least equal to array->n.
 */
void driver_affect_real(double *x, driver_array *array)
{
  int i;

  for (i = 0; i < array->n; i++) {
    switch (array->array[i].type) {
    case INTEGER:
      x[i] = array->array[i].d.i;
      break;
    case REAL:
      x[i] = array->array[i].d.d;
      break;
    case STRING:
      x[i] = array->array[i].d.s[0];
      break;
    }
  }
}


void driver_affect_integer(int *x, driver_array *array)
{
  int i;

  for (i = 0; i < array->n; i++) {
    switch (array->array[i].type) {
    case INTEGER:
      x[i] = array->array[i].d.i;
      break;
    case REAL:
      x[i] = (int)(array->array[i].d.d);
      break;
    case STRING:
      x[i] = array->array[i].d.s[0];
      break;
    }
  }
}

void *driver_array_address(AMLET_Driver *dr, int id, int pos, int *type)
{
  /* Could be better designed ? */
  *type = INTEGER;

  switch (id) {
  case ALTTYPES:
    *type = INTEGER;
    return &(dr->alt_types[pos]);
  case ALTVALUES:
    *type = REAL;
    return &(dr->alt_values[pos]);
  default:
    c_warning("invalid array");
    break;
  }

  return 0;
}

void *driver_get_address(int id)
{
  return NULL;
}

/**
 * Assign the read values to the model.
 *
 * @param dr is the built driver.
 * @param m is an initialized AMLET model.
 */
void driver_assign_model(AMLET_Driver *dr, AMLET_Model *m)
{
  OP_Problem *p = amlet_get_problem(m);

  if (dr->var[NAME].s)
    op_problem_set_name(p, dr->var[NAME].s);
  else
    op_problem_set_name(p, AMLET_DEFAULT_PROBLEM_NAME);
  m->sampling_id = dr->var[SAMPLING].i;
}

/**
 * Parse the driver file.
 *
 * @param f is the file to parse.
 * @param dr is a pointer towards a driver structure; it is
 *   uninitialized on input, and is on output fullfilled with the
 *   information that has been read in the file f, if this file has the
 *   correct format.
 *
 * @return 0 if the file has been correctly  parsed, another integer otherwise
 */
int amlet_parse_driver(FILE *f, AMLET_Driver *dr)
{
  int out;

  /* Initialize the driver */
  driver_init_defaults(dr);
  /* Parse the driver file */
  yyinit(f);
  amlet_parser_init(dr);
  out = yyparse();
  amlet_parser_clean();

  return out;
}

/**
 * Create a filled driver structure based on a driver file.
 *
 * @param name is the name of the driver file.
 *
 * @return NULL if the driver has not been correctly created, in case
 * the driver file cannot be open or it cannot be correctly parsed.
 */
AMLET_Driver *amlet_driver(const char *name)
{
  AMLET_Driver *dr = NULL;

  FILE *f = fopen(name, "r");
  if (f) {
    dr = _malloc(AMLET_Driver); 

    if (amlet_parse_driver(f, dr)) {
      driver_clean(dr);
      free(dr);
      dr = NULL;
    }
    fclose(f);
  }

  return dr;
}
