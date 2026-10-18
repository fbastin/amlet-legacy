/**
 * driver.h
 *
 * AMLET library
 * Copyright (C) 2002-2010 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#ifndef __AMLET_LIBRARY_DRIVER_H__
#define __AMLET_LIBRARY_DRIVER_H__

#include <oratio/ctools.h>
#include <oratio/moca.h>

#include <amlet/amlet_config.h>
#include <amlet/dcl.h>
#include <amlet/driver_parser.h>
#include <amlet/model.h>

#ifdef __cplusplus
extern "C" {
#endif

  /* array identifiers */

  /* scalar identifiers */
  enum {
    NP,              /* number of people */
    NALT,            /* number of alternatives */
    NOBS,            /* number of observations */
    NVAR,
    IDDEP,           /* number of the variable that is the dependent
		        variable */
    CENSOR,
    IDCENSOR,

    RESCALE,         /* 1 if some independent variables has to be rescale,
			0 otherwise. */
    BIAS,            /* 1 is we apply the bias correction, 0 else. */

    NITER,           /* number of iterations allowed in the maximization
                        process */
    EPS,             /* tolerance for convergence */
    SEED1,           /* seed for random generator */
    NREP,            /* maximum sampling size */
    VARREP,          /* variable sample size indicator */
    LEVEL,           /* signification level for statistical tests */
    LEVEL_ACCURACY,  /* signification level for objective function
			evaluation */

    FORMAT,          /* data format */

    METHOD,          /* algorithm used for the optimization process */
    HESSIAN_UPDATE,  /* update strategy for the Hessian update */
    HESSIAN_APPROX,  /* use approximation instead of the true Hessian when
			computing statistical informations */

    RANGEN,          /* random numbers generator identifier */
    LEGACY_RNG,      /* 1 to reproduce the MRG32k3a draws of ORATIO 0.7.5
			(AMLET 0.13.1), 0 for L'Ecuyer's MRG32k3a */
    SAMPLING,        /* sampling method */
    RANDOMIZATION,   /* randomization method */
    QMCREP,          /* number of samplings for accuracy estimation with
			quasi-Monte Carlo draws */

    VALIDATION,      /* validation indicator */
    STATISTICS,      /* compute error and bias of simulation


		        (automatically done if MCBTRDA method is selected) */

    ROBUST,          /* robust estimation of the information matrix */

    /* sampling */
    DRAWS,
    ALTERNATIVES_SAMPLING,

    /* data information */
    SYNTHETIC,       /* true if we use synthetic data, false else */

    /* Nests */
    NNESTS,          /* number of nests */

    /* strings */
    NAME,

    /* filenames */
    DRAWSNAME,

    END_RESERVED_WORDS
  };

  /* Random and quasi-random numbers generators */

  static const c_symbol randomization_halton[] = {
    {"none", MOCA_HALTON_RANDOMIZATION_NONE},
    {"scrambling", MOCA_HALTON_RANDOMIZATION_SCRAMBLING},
    {"shifting", MOCA_HALTON_RANDOMIZATION_SHIFTING},
    {"shifting_scrambling", MOCA_HALTON_RANDOMIZATION_SHIFTING_SCRAMBLING},
    {"shuffling", MOCA_HALTON_RANDOMIZATION_SHUFFLING},
    {"shifting_shuffling", MOCA_HALTON_RANDOMIZATION_SHIFTING_SHUFFLING}
  };

  enum {
    BTYPE = 0,     /* type of parameters to estimate */
    BGENTYPE,
    B,
    BGEN,
    IGNORE,
    RESCALEIDX,
    RESCALEVAL,
    END_ARRAY_WORDS
  };

  enum {
    NESTS = 0,
    END_OTHERS
  };

  /* reserved words */
  static const c_symbol reserved_words[] = {
    {"np", NP},
    {"nalt", NALT},
    {"nobs", NOBS},
    {"nvar", NVAR},
    {"iddep", IDDEP},

    {"censor", CENSOR},
    {"idcensor", IDCENSOR},

    {"rescale", RESCALE},

    {"bias", BIAS},

    {"niter", NITER},
    {"eps", EPS},

    /* Random numbers and draws generator parameters */
    {"seed1", SEED1},
    {"nrep", NREP},
    {"varrep", VARREP},
    {"level", LEVEL},
    {"level_accuracy", LEVEL_ACCURACY},

    {"format", FORMAT},

    {"method", METHOD},
    {"hessian_update", HESSIAN_UPDATE},
    {"hessian_approx", HESSIAN_APPROX},

    {"rangen", RANGEN},
    {"legacy_rng", LEGACY_RNG},
    {"sampling", SAMPLING},
    {"randomization", RANDOMIZATION},
    {"qmcrep", QMCREP},

    {"validation", VALIDATION},
    {"robust", ROBUST},
    {"statistics", STATISTICS},

    {"draws", DRAWS},
    {"alternatives_sampling", ALTERNATIVES_SAMPLING},

    {"synthetic", SYNTHETIC},

    {"nnests", NNESTS},

    {"problem", NAME},
    {"drawsname", DRAWSNAME},

    /* arrays */
    {"btype", BTYPE},
    {"bgentype", BGENTYPE},
    {"b", B},
    {"bgen", BGEN},
    {"ignore", IGNORE},
    {"rescalidx", RESCALEIDX},
    {"rescalval", RESCALEVAL},

    /* others */ 
    {"nests", NESTS}
  };

  /* double entry arrays */

  enum {
    ALTTYPES = 0,
    ALTVALUES,
    END_DBL_ARRAY_WORDS
  };

  /* arrays */
  static const c_symbol arrays[] = {
    {"alttypes", ALTTYPES},
    {"altvalues", ALTVALUES},
  };

  /* optimization methods */

  /* In the order of the enum in dcl.h: io.c prints methods[method].name
     (AMLET <= 0.13.1 printed mcrtrda for mcbtrda and vice versa). */
  static const c_symbol methods[] = {
    {"none", NO_OPTIMIZATION},
    {"mcbtr", MCBTR},
    {"mcrtr", MCRTR},
    {"mcrtrda", MCRTRDA},
    {"mcbtrda", MCBTRDA},
    {"mcbfgs", MCBFGS},
    {"mcbfgsda", MCBFGSDA},
    {"mcbtrc", MCBTRC}
  };

  static const c_symbol hessian_updates[] = {
    {"bfgs", HESSIAN_BFGS},
    {"bhhh", HESSIAN_BHHH},
    {"sr1", HESSIAN_SR1}
  };

  static const c_symbol rangen[] = {
    {"ecuyer", RAN_LECUYER88B},
    {"mrg32k3a", RAN_MRG32K3A},
    {"standard", RAN_STANDARDMINIMAL},
    {"parkmiller", RAN_STANDARDMINIMAL},
    {"mt19937", RAN_MT19937}
  };

  typedef union {
    int    i;
    double d;
    char  *s;
  } driver_data;

  typedef struct {
    driver_data d;
    int type;
  } driver_var;

  typedef struct {
    driver_var *array;
    int n;
  } driver_array;

#define RESERVED_WORDS NP
#define ARRAY_WORDS    B

typedef struct {
    driver_data var[END_RESERVED_WORDS];
    driver_array array[END_ARRAY_WORDS];

    /* TODO: put them in the reserved words */
    int generation;
    char *obsname; /* name of the observations file */
    char *output;  /* name of the output file */

    driver_array *nests;   /* nests gestion */

    /* for the generation */
    int **alt_types;
    double **alt_values;
    int *alt_ntypes;    /* number of values given in alttypes[i] */
    int *alt_nvalues;   /* number of values given in altvalues[i] */
  } AMLET_Driver;

  /**
   * Assign the read values to the model.
   *
   * @param dr is the built driver.
   * @param m is an initialized AMLET model.
   */
  void driver_assign_model(AMLET_Driver *dr, AMLET_Model *m);


  /* Lexical analyzer functions */
  void yyinit(FILE *f);

  void driver_init_defaults(AMLET_Driver *dr);

  void driver_post_parsing(AMLET_Driver *dr, dcl_data *md);

  void driver_clean(AMLET_Driver *dr);

  /**
   * Check the parameters defined in the driver
   *
   * @param dr is a pointer to the driver structure.
   * @param data is a pointer to the structure containing important data about
   *   the problem under consideration.
   */
  void driver_check_parameters(AMLET_Driver *dr, amlet_config *ac,
			       dcl_data *data);

  int driver_check_array_pos(AMLET_Driver *dr, int id, int pos);

  /**
   * Do some checks after the parsing
   */
  void driver_check(AMLET_Driver *dr);

  /***************************************************************************/

  void driver_set_generation_process(AMLET_Driver *dr, char *obsname,
				     int nalt);

  /***************************************************************************/

  /**
   * @param x is a double vector of size at least equal to array->n.
   */
  void driver_affect_real(double *x, driver_array *array);
  void driver_affect_integer(int *x, driver_array *array);

  void driver_free_array(AMLET_Driver *dr, int id);

  void *driver_array_address(AMLET_Driver *dr, int id, int pos, int *type,
			     int **count);

  /**
   * Assign a value read in the driver file to a reserved word.
   *
   * @param word is the position of the word in reserved_words[].
   * @param type is the type of the value read: INTEGER, REAL or STRING.
   */
  void driver_assign_value(AMLET_Driver *dr, int word, int type,
			   driver_data value);

  /**
   * Assign a list of values read in the driver file to a reserved word.
   *
   * @param word is the position of the word in reserved_words[].
   */
  void driver_assign_array(AMLET_Driver *dr, int word,
			   const driver_var *values, int n);

  void driver_array_action(AMLET_Driver *dr, void *address, int n);

  void *driver_get_address(int id);

  /***************************************************************************/

  void driver_identifier_action(AMLET_Driver *dr, int id);

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
  int amlet_parse_driver(FILE *f, AMLET_Driver *dr);

  /**
   * Create a filled driver structure based on a driver file.
   *
   * @param name is the name of the driver file.
   *
   * @return NULL if the driver has not been correctly created, in case
   * the driver file cannot be open or it cannot be correctly parsed.
   */
  AMLET_Driver *amlet_driver(const char *name);

#ifdef __cplusplus
}
#endif

#endif  /* __AMLET_LIBRARY_DRIVER_H__ */
