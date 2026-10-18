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
#include "lattice.h"

#include "model.h"

#include <oratio/ctools.h>

#include <limits.h>
#include <math.h>
#include <string.h>

#define DRIVER_DEFAULT_NREP 1000

/**
 * External parsing functions, defined in the Yacc source code driver_parser.y
 */

extern void amlet_parser_init(AMLET_Driver *d);
extern void amlet_parser_clean();
extern int yyparse ();
extern int lineno;

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
  dr->alt_ntypes = 0;
  dr->alt_nvalues = 0;

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
    free(dr->alt_ntypes);
    free(dr->alt_nvalues);
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

  /* The estimation code has no nests (exec.c sets data.n to NULL); with
     nnests > 1, AMLET <= 0.13.1 wrote past its arrays. */
  if (dr->var[NNESTS].i != 1)
    c_error("nested logit models are not supported: nnests must be 1");

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

    /* mcbtrc constrains the spline coefficients to be monotonic; without
       any, AMLET <= 0.13.1 crashed. The method is then mcbtr. */
    if (dr->var[METHOD].i == MCBTRC && ds_list_length(data->np_list) == 0) {
      amlet_warning(stderr, 1, "mcbtrc constrains spline coefficients, and the model has none; mcbtr is used instead");
      dr->var[METHOD].i = MCBTR;
    }

    /* Faure sequences are not implemented in the draws generation:
       AMLET <= 0.13.1 used Monte Carlo draws instead. */
    if (dr->var[SAMPLING].i == MOCA_FAURE)
      c_error("Faure sequences are not available; use sampling = halton, sobol, lattice, mlhs or monte_carlo");

    /* Rank-1 lattices: use the number of points of a tabulated lattice
       (a prime number), not more than nrep. AMLET <= 0.13.1 kept nrep
       points with the generating vector of another lattice. */
    if (dr->var[SAMPLING].i == MOCA_RANK1_LATTICE) {
      int *l = amlet_get_default_lattice();
      int inc = amlet_get_default_lattice_max_dimension()+1;
      int k, size = l[0];

      if (data->n_random_coeff > amlet_get_default_lattice_max_dimension())
	amlet_warning(stderr, 1, "rank-1 lattices are available up to dimension %d; Monte Carlo draws are used",
		      amlet_get_default_lattice_max_dimension());
      else {
	for (k = 0; k < amlet_get_default_lattice_number(); k++)
	  if (l[k*inc] <= dr->var[NREP].i) size = l[k*inc];
	if (size != dr->var[NREP].i) {
	  amlet_warning(stderr, 1, "rank-1 lattice: nrep set to %d, the number of points of the closest tabulated lattice", size);
	  dr->var[NREP].i = size;
	}
      }
    }

    /* The variable sample size methods estimate the simulation error as
       for independent draws: with quasi-Monte Carlo draws, they stall or
       stop with too few draws. Use the method with a fixed sample size. */
    if (dr->var[SAMPLING].i != MOCA_MONTE_CARLO
	&& (dr->var[METHOD].i == MCBTRDA || dr->var[METHOD].i == MCRTRDA)) {
      int fixed = (dr->var[METHOD].i == MCBTRDA) ? MCBTR : MCRTR;

      amlet_warning(stderr, 1, "%s varies the number of draws, which requires Monte Carlo draws (sampling = monte_carlo); %s is used instead, with %d draws",
		    methods[dr->var[METHOD].i].name, methods[fixed].name,
		    dr->var[NREP].i);
      dr->var[METHOD].i = fixed;
    }

    /* The correction of the simulation bias also estimates the error
       of the simulated probabilities as for independent draws. */
    if (dr->var[SAMPLING].i != MOCA_MONTE_CARLO && dr->var[BIAS].i) {
      amlet_warning(stderr, 1, "the correction of the simulation bias requires Monte Carlo draws (sampling = monte_carlo); bias is set to 0");
      dr->var[BIAS].i = data->bias = 0;
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
	int k, nvalues = 0;

	if (!dr->alt_types[i] || !dr->alt_values[i]) {
	  c_error("missing definition for alternative %d\n", i);
	}
	/* AMLET <= 0.13.1 read past shorter lists. */
	if (dr->alt_ntypes[i] < data->nvar)
	  c_error("alttypes[%d] has %d values; the model has %d variables",
		  i, dr->alt_ntypes[i], data->nvar);
	for (k = 0; k < data->nvar; k++) {
	  if (amlet_is_nonparametric(dr->alt_types[i][k]))
	    c_error("alttypes[%d]: nonparametric distributions cannot be used to generate the variables", i);
	  nvalues += moca_count_parameters(dr->alt_types[i][k], NULL);
	}
	if (dr->alt_nvalues[i] < nvalues)
	  c_error("altvalues[%d] has %d values; alttypes[%d] needs %d",
		  i, dr->alt_nvalues[i], i, nvalues);
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
  /* AMLET <= 0.13.1 crashed, or divided by zero */
  if (dr->var[NALT].i < 2)
    c_error("nalt (the number of alternatives) must be given, and at least 2");
  if (dr->var[NOBS].i < 1)
    c_error("nobs (the number of observations) must be at least 1");
  if (dr->var[NVAR].i - (dr->var[CENSOR].i ? 1 : 0)
      - (dr->var[IDDEP].i ? 1 : 0) < 1)
    c_error("nvar = %d leaves no explanatory variable", dr->var[NVAR].i);

  if (!dr->var[CENSOR].i) {
    if (dr->var[IDCENSOR].i) {
      c_warning("IDCENSOR set while CENSOR equal to 0; IDCENSOR reset to 0\n");
      dr->var[IDCENSOR].i = 0;
    }
  }
  /* AMLET <= 0.13.1 then left the availabilities undefined */
  else if (dr->var[IDCENSOR].i < 1 || dr->var[IDCENSOR].i > dr->var[NVAR].i)
    c_error("censor = 1: idcensor must give the row of the availabilities, from 1 to nvar = %d",
	    dr->var[NVAR].i);
  if (dr->var[NITER].i < 0) {
    c_warning("number of iterations set to a negative value; reset to default value of %d\n", AMLET_DEFAULT_NITER);
    dr->var[NITER].i = AMLET_DEFAULT_NITER;
  }
  /* AMLET <= 0.13.1 seeded the generator with 2^64 - |seed1|, not a valid
     MRG32k3a state */
  if (dr->var[SEED1].i < 0)
    c_error("seed1 must be positive (or 0, for a seed taken from the clock)");
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

  if (nalt <= 0)
    c_error("line %d: nalt must be given before generate", lineno);

  if (dr->alt_types) {
    c_warning("generation process redefined\n");
    for (i = 0; i < nalt; i++) free(dr->alt_types[i]);
    for (i = 0; i < nalt; i++) free(dr->alt_values[i]);
    free(dr->alt_types);
    free(dr->alt_values);
    free(dr->alt_ntypes);
    free(dr->alt_nvalues);
  }

  dr->obsname = obsname;
  /* AMLET <= 0.13.1 left them undefined: a missing alternative was
     detected only when malloc() returned zeroed memory. */
  dr->alt_types = _mallocs(nalt, int *);
  dr->alt_values = _mallocs(nalt, double *);
  dr->alt_ntypes = _mallocs(nalt, int);
  dr->alt_nvalues = _mallocs(nalt, int);
  for (i = 0; i < nalt; i++) {
    dr->alt_types[i] = NULL;
    dr->alt_values[i] = NULL;
    dr->alt_ntypes[i] = dr->alt_nvalues[i] = 0;
  }
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

void *driver_array_address(AMLET_Driver *dr, int id, int pos, int *type,
			   int **count)
{
  /* Could be better designed ? */
  *type = INTEGER;

  switch (id) {
  case ALTTYPES:
    *type = INTEGER;
    *count = &(dr->alt_ntypes[pos]);
    return &(dr->alt_types[pos]);
  case ALTVALUES:
    *type = REAL;
    *count = &(dr->alt_nvalues[pos]);
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
 * Type of the value of a scalar reserved word: INTEGER, REAL or STRING.
 */
static int driver_value_type(int id)
{
  switch (id) {
  case EPS:
  case LEVEL:
  case LEVEL_ACCURACY:
    return REAL;
  case NAME:
  case DRAWSNAME:
    return STRING;
  default:
    return INTEGER;
  }
}

/**
 * Assign a value read in the driver file to a reserved word.
 *
 * AMLET <= 0.13.1 stored the value in the field of its own type, whatever
 * the type of the parameter: niter = 100.0 gave 0 iterations, eps = 1 left
 * eps almost unchanged, and nrep = {1000} or b = 0 wrote past the
 * parameter.
 *
 * @param word is the position of the word in reserved_words[].
 * @param type is the type of the value read: INTEGER, REAL or STRING.
 */
void driver_assign_value(AMLET_Driver *dr, int word, int type,
			 driver_data value)
{
  int nwords = sizeof(reserved_words)/sizeof(reserved_words[0]);
  int id = reserved_words[word].token;
  char *name = reserved_words[word].name;

  if (word >= nwords - END_OTHERS) {
    /* nests = n is nnests = n in AMLET <= 0.10; it crashed AMLET 0.13.1 */
    c_warning("line %d: nests is the name of nnests in AMLET <= 0.10; use nnests",
	      lineno);
    id = NNESTS;
    name = "nnests";
  }
  else if (word >= nwords - END_ARRAY_WORDS - END_OTHERS)
    c_error("line %d: %s expects a list of values, as in %s = {...}",
	    lineno, name, name);

  switch (driver_value_type(id)) {
  case REAL:
    if (type == STRING)
      c_error("line %d: %s expects a number", lineno, name);
    dr->var[id].d = (type == INTEGER) ? value.i : value.d;
    break;
  case STRING:
    if (type != STRING)
      c_error("line %d: %s expects a string, as in %s = \"...\"",
	      lineno, name, name);
    free(dr->var[id].s);
    dr->var[id].s = value.s;
    break;
  default:
    if (type == REAL && value.d == rint(value.d) && fabs(value.d) <= INT_MAX)
      value.i = (int)value.d;
    else if (type != INTEGER)
      c_error("line %d: %s expects an integer", lineno, name);
    /* the arrays of the generation process have nalt elements */
    if (id == NALT && dr->alt_types && value.i != dr->var[NALT].i)
      c_error("line %d: nalt must be given before generate", lineno);
    /* parsed, but no code uses them */
    if ((id == RESCALE || id == VARREP || id == SYNTHETIC) && value.i) {
      c_warning("line %d: %s is not implemented in AMLET 0.13.1; it is ignored",
		lineno, name);
      value.i = 0;
    }
    dr->var[id].i = value.i;
  }
}

/**
 * Assign a list of values read in the driver file to a reserved word.
 *
 * @param word is the position of the word in reserved_words[].
 */
void driver_assign_array(AMLET_Driver *dr, int word,
			 const driver_var *values, int n)
{
  int nwords = sizeof(reserved_words)/sizeof(reserved_words[0]);
  int id = reserved_words[word].token;
  driver_array *a = dr->array + id;

  if (word < nwords - END_ARRAY_WORDS - END_OTHERS
      || word >= nwords - END_OTHERS)
    c_error("line %d: %s expects a single value, not a list",
	    lineno, reserved_words[word].name);
  /* parsed, but no code uses them */
  if ((id == IGNORE || id == RESCALEIDX || id == RESCALEVAL) && n)
    c_warning("line %d: %s is not implemented in AMLET 0.13.1; it is ignored",
	      lineno, reserved_words[word].name);

  free(a->array);
  a->n = n;
  a->array = (!n) ? NULL : _mallocs(n, driver_var);
  if (n) memcpy(a->array, values, n*sizeof(driver_var));
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
