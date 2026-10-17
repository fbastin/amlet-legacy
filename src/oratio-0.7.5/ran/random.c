/**
 * random.c
 *
 * ORATIO Library
 * Copyright (C) 1999-2007 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include "oratio/ran.h"

#include <time.h>
#include <string.h>

#include "lecuyer88a.h"
#include "lecuyer88b.h"
#include "mrg32k3a.h"
#include "mt19937.h"
#include "standardminimal.h"

#include <oratio/oratio.h>

static const c_symbol generator_name[] = {
  { "LEcuyer88a", RAN_LECUYER88A},
  { "LEcuyer88b", RAN_LECUYER88B},
  { "MRG32K3A", RAN_MRG32K3A},
  { "MT19937", RAN_MT19937},
  { "Standard_Minimal", RAN_STANDARDMINIMAL}
};

/**
 * Store a copy of the seed vector inside the random numbers generator.
 *
 * @param random is an initialized random numbers generator, except
 *   that the seed vector is undefined on input.
 * @param seed is the initial seed vector that has to be stored in
 *   random.
 * @param n is the dimension of seed.
 */
void ran_random_set_seed(Random *random, unsigned long *seed, int n)
{
  if (seed) {
    random->seed = _reallocs(random->seed, n, unsigned long);
    memcpy(random->seed, seed, n*sizeof(unsigned long));
  }
  else
    random->seed = NULL;
}

/**
 * Create a new generator with a initial seed.
 *
 * @param seed is the seed vector used to initialize the random
 *   numbers generator. A copy of this vector is stored inside the
 *  generator.
 * @param n is the dimension of seed.
 *
 * @return the newly created random object.
 */
Random *ran_random_new(int id, unsigned long *seed, int n)
{
  Random *random = _malloc(Random);

  random->seed = random->rng = random->parameters = NULL;
  random->id = (id) ? id : RAN_DEFAULT;
  random->type = RAN_UNIFORM;

  random->l = 0.0;
  random->h = 1.0;
  random->p = 0.0;

  ran_random_set_generator(random, random->id, seed, n);

  return random;
}

/**
 * Create a new generator using default values.
 *
 * @return the newly created random object.
 */
Random *ran_random()
{
  return ran_random_new(RAN_DEFAULT, NULL, 0);
}

/**
 * Free the memory associated to a the underlying U(0,1)-random numbers
 * generator.
 *
 * @param random is a previously initialized random generator on
 *   input, and is undefined on output, the associated memory having
 *   been freed.
 */
void ran_random_free_rng(Random *random)
{
  if (random->rng) {
    switch (random->id) {
    case RAN_LECUYER88A:
      ran_lecuyer88a_free(random->rng);
      break;
    case RAN_LECUYER88B:
      ran_lecuyer88b_free(random->rng);
      break;
    case RAN_MRG32K3A:
      ran_mrg32k3a_free((random->rng));
      break;
    case RAN_MT19937:
      ran_mt19937_free(random->rng);
      break;
    case RAN_STANDARDMINIMAL:
      ran_standardminimal_free(random->rng);
      break;
    default:
      break;
    }
    random->rng = NULL;
  }
}

/**
 * Free the memory associated to a random numbers generator.
 */
void ran_random_free(Random *random)
{
  if (random) {
    ran_random_free_rng(random);
    if (random->seed) free(random->seed);
    if (random->parameters) free(random->parameters);
    free(random);
  }
}

/**
 * Uniform Generator Management
 */

/**
 * Set the uniform random numbers generator associated to a random object.
 *
 * @param random is an initialized random numbers generator.
 * @param id is the generator identifier.
 * @param seed is the seed vector used to initialize the random
 *   numbers generator.
 * @param n is the dimension of seed.
 */
void ran_random_set_generator(Random *random, const int id,
			      unsigned long *seed, int n)
{
  ran_random_free_rng(random);
  if (n) ran_random_set_seed(random, seed, n);

  switch (id) {
  case RAN_LECUYER88A:
    random->rng = ran_lecuyer88a_new(seed, n);
    random->id = RAN_LECUYER88A;
    random->rng_get_val = ran_lecuyer88a_get_val;
    break;
  case RAN_LECUYER88B:
    random->rng = ran_lecuyer88b_new(seed, n);
    random->id = RAN_LECUYER88B;
    random->rng_get_val = ran_lecuyer88b_get_val;
    break;
    /*  case RAN_ECUYERSIMARD99:
    random->rng = ran_ecuyersimard99_new(seed, n);
    random->id = RAN_ECUYERSIMARD99;
    random->rng_get_val = ran_ecuyersimard99_get_val;
    breal;
    */
  case RAN_MRG32K3A:
    random->rng = ran_mrg32k3a_new(seed, n);
    random->id = RAN_MRG32K3A;
    random->rng_get_val = ran_mrg32k3a_get_val;
    break;
  case RAN_MT19937:
    random->rng = ran_mt19937_new(seed, n);
    random->id = RAN_MT19937;
    random->rng_get_val = ran_mt19937_get_val;
    break;
  case RAN_STANDARDMINIMAL: 
  default: 
    random->rng =
      (seed) ? ran_standardminimal_new(*seed) :  ran_standardminimal_new(0);
    random->id = RAN_STANDARDMINIMAL;
    random->rng_get_val = ran_standardminimal_get_val;
    break;
  }
}

/**
 * Return the generator name.
 *
 * @param random is an initialized random numbers generator.
 *
 * @return a pointer to the name of the generator. There is no new
 *   memory allocation, so the returned pointer has not to be freed.
 */
char *ran_random_generator_name(Random *random)
{
  return generator_name[random->id-1].name;
}

/**
 * Return the number of parameters needed for a generator of type 'type'
 *
 * @param random is an initialized random numbers generator.
 *
 * @return the number of parameters used to described the distribution
 */
int ran_random_get_parameters_number(Random *random)
{
  switch (random->type) {
  case RAN_ARCSINE:
    return RAN_ARCSINE_NPAR;
  case RAN_BERNOUILLI:
    return RAN_BERNOUILLI_NPAR;
  case RAN_CAUCHY:
    return RAN_CAUCHY_NPAR;
  case RAN_EXPONENTIAL:
    return RAN_EXPONENTIAL_NPAR;
  case RAN_GAMMA:
    return RAN_GAMMA_NPAR;
  case RAN_GUMBEL:
    return RAN_GUMBEL_NPAR;
  case RAN_GEOMETRIC:
    return RAN_GEOMETRIC_NPAR;
  case RAN_INTEGER:
    return RAN_INTEGER_NPAR;
  case RAN_JOHNSONSB:
    return RAN_JOHNSONSB_NPAR;
  case RAN_JOHNSONSU:
    return RAN_JOHNSONSU_NPAR;
  case RAN_LAPLACE:
    return RAN_LAPLACE_NPAR;
  case RAN_LOGNORMAL:
    return RAN_LOGNORMAL_NPAR;
  case RAN_NORMAL:
    return RAN_NORMAL_NPAR;
  case RAN_STUDENT:
    return RAN_STUDENT_NPAR;
  case RAN_UNIFORM:
    return RAN_UNIFORM_NPAR;
  case RAN_TRIANGULAR:
    return RAN_TRIANGULAR_NPAR;
  case RAN_TRUNCATED_NORMAL:
    return RAN_TRUNCATED_NORMAL_NPAR;
  default:
    return 0;
  };
}

/**
 * State gestion.
 */

/**
 * Return the current state of the generator.
 *
 * @param random is an initialized random numbers generator.
 *
 * @return the state associated to the random numbers generators, or
 *   NULL if the generator type is unidentified.
 */
void *ran_random_get_state(Random *random)
{
  void *state = NULL;

  switch (random->id) {
  case RAN_LECUYER88A:
    state = ran_lecuyer88a_get_state(random->rng);
    break;
  case RAN_LECUYER88B:
    state = ran_lecuyer88b_get_state(random->rng);
    break;
  case RAN_MRG32K3A:
    state = ran_mrg32k3a_get_state(random->rng);
    break;
  case RAN_MT19937:
    state = ran_mt19937_get_state(random->rng);
    break;
  case RAN_STANDARDMINIMAL:
    state = ran_standardminimal_get_state(random->rng);
    break;
  default:
    break;
  }

  return state;
}

/**
 * Copy the generator current state.
 *
 * @param random is an initialized random numbers generator.
 * @param n contains on output the state length.
 *
 * @return a newly allocated pointer to a copy of the state of the
 *   generator; the state is store as a unsigned long vector.
 */
unsigned long *ran_random_copy_state(Random *random, int *n)
{
  void *state = ran_random_get_state(random);
  unsigned long *s = NULL;

  int i;

  switch (random->id) {
  case RAN_LECUYER88A:
  case RAN_LECUYER88B:
    *n = 2;
    break;
  case RAN_MRG32K3A:
    *n = 6;
    break;
  case RAN_MT19937:
    *n = MT19937_STATE_LENGTH;
    break;
  case RAN_STANDARDMINIMAL:
  default:
    *n = 1;
    break;
  }

  s = _mallocs(*n, unsigned long);
  if (random->id == RAN_MRG32K3A)
    for (i = 0; i < 6; i++) s[i] = ((double *)state)[i];
  else
    memcpy(s, state, *n*sizeof(unsigned long));

  return s;
}

/**
 * Set the current state of the generator.
 *
 * @param random is an initialized random numbers generator.
 * @param state is the new state to associate to the random numbers
 *   generator. The precise type must correspond to the associated
 *   random object. If state is NULL, nothing if performed.
 */
void ran_random_set_state(Random *random, void *state)
{
  switch (random->id) {
  case RAN_LECUYER88A:
    ran_lecuyer88a_set_state((LEcuyer88a *)random, state);
    break;
  case RAN_LECUYER88B:
    ran_lecuyer88b_set_state((LEcuyer88b *)random, state);
    break;
  case RAN_STANDARDMINIMAL:
    ran_standardminimal_set_state(random->rng, state);
    break;
  case RAN_MRG32K3A:
    ran_mrg32k3a_set_state(random->rng, state);
    break;
  case RAN_MT19937:
    ran_mt19937_set_state(random->rng, state);
    break;
  default:
    c_warning("unrecognized method identifier for random numbers generation");
    break;
  }
}

/**
 * Create a new generator with an initial seed and the parameters
 * vector size.
 *
 * @param seed is the seed vector used to initialize the random
 *   numbers generator. A copy of this vector is stored inside the
 *  generator.
 * @param n is the dimension of seed.
 * @param size is the size in bytes of the parameters vector guiding
 *   the generator behavior.
 *
 * @return the newly created random object.
 */
Random *ran_random_new_with_size(unsigned long *seed, int n,
				 size_t size)
{
  Random *random = ran_random_new(0, seed, n);

  random->parameters = (size) ? malloc(size) : NULL;

  return random;
}

/**
 * Draw an uniform random number in (F(a), F(b)).
 *
 * @param random is an initialized random numbers generator.
 */
double ran_random_get_val(const Random *random)
{
  double low = random->l-random->p;
  return low+random->rng_get_val(random->rng)*(random->h-low);
}

/**
 * Return a draw from an univariate random number generator, using the type
 * associated to the generator.
 *
 * @param random is an initialized random numbers generator.
 *
 * @return a draw from the univariate random number generator.
 */
double ran_random_next(const Random *random)
{
  switch (random->type) {
  case RAN_ARCSINE:
    return ran_arcsine_get_val(random);
  case RAN_BERNOUILLI:
    return ran_bernouilli_get_val(random);
  case RAN_CAUCHY:
    return ran_cauchy_get_val(random);
  case RAN_DIRAC:
    return ran_dirac_get_val(random);
  case RAN_EXPONENTIAL:
    return ran_exponential_get_val(random);
  case RAN_GAMMA:
    return ran_gamma_get_val(random);
  case RAN_GEOMETRIC:
    return ran_geometric_get_val(random);
  case RAN_GUMBEL:
    return ran_gumbel_get_val(random);
  case RAN_INTEGER:
    return ran_intrandom_get_val(random);
  case RAN_JOHNSONSB:
    return ran_johnsonSB_get_val(random);
  case RAN_JOHNSONSU:
    return ran_johnsonSU_get_val(random);
  case RAN_LAPLACE:
    return ran_laplace_get_val(random);
  case RAN_LOGNORMAL:
    return ran_lognormal_get_val(random);
  case RAN_NORMAL:
    return ran_normal_get_val(random);
  case RAN_SP_SPLINE:
    return ran_sp_spline_get_val(random);
  case RAN_STUDENT:
    return ran_student_get_val(random);
  case RAN_TRIANGULAR:
    return ran_triangular_get_val(random);
  case RAN_TRUNCATED_NORMAL:
    return ran_truncatednormal_get_val(random);
  case RAN_UNIFORM:
    return ran_random_get_val(random);
  default:
    return 0;
  };
}

void ran_random_set_parameters(Random *random,
			       void *parameters, size_t size)
{
  if (random) {
    if (random->parameters) free(random->parameters);
    if (size) {
      random->parameters = malloc(size);
      memcpy(random->parameters, parameters, size);
    }
    else
      random->parameters = NULL;
  }
}

/**
 * Define the parameters allowing truncation when inversion
 * technique is used.
 *
 * @param random is an initialized random numbers generator.
 * @param l is a lower bound on the untruncated cdf.
 * @param h is an upper bound on the untruncated cdf.
 * @param p is the value of P[X=a], where a is a real such that l =
 *   F(a).
 */
void ran_random_set_truncation(Random *random, double l,
				      double h, double p)
{
  random->l = l;
  random->h = h;
  random->p = p;
}

/**
 * Produce an integer, using the intern clock, in (0,2^31-1).
 * This function is very poor...
 */
unsigned long ran_init_random()
{
  time_t    timer;
  struct tm *tblock;
  long      x, y, z;

  long im1 = 2147483563;       /* 2^31-1 */

  timer = time(NULL);          /* get the hour */
  tblock = localtime(&timer);  /* convert the hour in structure */

  x = ((long)floor((tblock->tm_mday*tblock->tm_min+1)/(tblock->tm_mon+1))+1)*
    (tblock->tm_year+1900)+tblock->tm_sec*(tblock->tm_min+1);

  /* This line ensure the the Schrage algorithm can be used. */
  if (x > sqrt(im1)) x = (long)floor(sqrt(x));

  y = nt_max(tblock->tm_sec+1, (tblock->tm_hour+tblock->tm_min+tblock->tm_sec)*
	     (tblock->tm_wday+tblock->tm_sec+1));

  z = nt_modular_schrage(x,y,im1);
  while (!z && y < im1-1) z = nt_modular_schrage(x,y++,im1);
  if (!z) z = 1;

  return x;
}

/**
 * Reset the current stream to its start.
 *
 * @param mrg is an initialized random numbers generator whose
 *   implementaion supports streams.
 */
void ran_random_resetstartstream(Random *random)
{
  if (random->rng) {
    switch (random->id) {
    case RAN_MRG32K3A:
      ran_mrg32k3a_resetstartstream((MRG32K3a *)(random->rng));
      break;
    default:
      fprintf(stderr, "Streams not supported by used generator!");
      break;
    }
  }
}

/**
 * Jump to the start of the next substream.
 *
 * @param mrg is an initialized random numbers generator whose
 *   implementaion supports streams.
 */
void ran_random_resetnextsubstream(Random *random)
{
  if (random->rng) {
    switch (random->id) {
    case RAN_MRG32K3A:
      ran_mrg32k3a_resetnextsubstream((MRG32K3a *)(random->rng));
      break;
    default:
      fprintf(stderr, "Streams not supported by used generator!");
      break;
    }
  }
}

/**
 * Jump to the state to the start of the current substream.
 *
 * @param mrg is an initialized random numbers generator whose
 *   implementaion supports streams.
 */
void ran_random_resetstartsubstream(Random *random)
{
  if (random->rng) {
    switch (random->id) {
    case RAN_MRG32K3A:
      ran_mrg32k3a_resetstartsubstream((MRG32K3a *)(random->rng));
      break;
    default:
      fprintf(stderr, "Streams not supported by used generator!");
      break;
    }
  }
}
