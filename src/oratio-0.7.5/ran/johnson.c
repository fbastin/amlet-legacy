/**
 * Johnson.c
 *
 * Random Numbers Library
 * Copyright (C) 2007 Cinzia Cirillo <ccirillo@umd.edu>,
 *                    Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

/**
 * Johnson variate.
 *
 * Reference: 
 * Slides from course "Simulation: aspects stochastiques" (IFT6561),
 * pp. 45-46
 * Pierre L'Ecuyer, Fabian Bastin, University of Montreal
 */

#include "oratio/ran.h"

/**
 * Create a new Johnson random variable generator with some initial seed.
 *
 * @param seed is the seed vector used to initialize the random
 *   numbers generator.
 * @param n is the dimension of the seed vector.
 * @param xi is the localisation parameter.
 * @param lamdda is the scale parameter.
 * @param gamma and delta are the form parameters.
 *
 * @return the newly allocated generator object.
 */
Random *ran_johnson_new(unsigned long *seed, int n,
			double xi, double lambda, double gamma, double delta)
{
  Random *random = ran_random_new_with_size(seed, n, 4*sizeof(double));

  ((double *)(random->parameters))[0] = xi;
  ((double *)(random->parameters))[1] = lambda;
  ((double *)(random->parameters))[2] = gamma;
  ((double *)(random->parameters))[3] = delta;

  return random; 
}

/**
 * Retrieve the parameters of a Johnson variate.
 *
 * @param xi contains on output the localisation parameter.
 * @param lamdda contains on ouput the scale parameter.
 * @param gamma and delta contain on output the form parameters.
 */
void ran_johnson_get_parameters(const Random *random,
				double *xi, double *lambda,
				double *gamma, double *delta)
{
  *xi     = ((double *)(random->parameters))[0];
  *lambda = ((double *)(random->parameters))[1];
  *gamma  = ((double *)(random->parameters))[2];
  *delta  = ((double *)(random->parameters))[3];
}

/**
 * Set the parameters of a Johnson random variate,
 *
 * @param xi is the localisation parameter.
 * @param lamdda is the scale parameter.
 * @param gamma and delta are the form parameters.
 */
void ran_johnson_set_parameters(Random *random,
				double xi, double lambda, 
				double gamma, double delta)
{  
  double temp[4];

  temp[0] = xi;
  temp[1] = lambda;
  temp[2] = gamma;
  temp[3] = delta;
  
  ran_random_set_parameters(random, temp, 4*sizeof(double));
}

/**
 * Unbounded Johnson.
 */

/**
 * Create a new unbounded Johnson random variable generator with some
 * initial state.
 *
 * @param seed is the seed vector used to initialize the random
 *   numbers generator.
 * @param n is the dimension of the seed vector.
 * @param xi is the localisation parameter.
 * @param lamdda is the scale parameter.
 * @param gamma and delta are the form parameters.
 *
 * @return the newly allocated generator object.
 */
Random *ran_johnsonSU_new(unsigned long *seed, int n,
				 double xi, double lambda,
				 double gamma, double delta)
{
  Random *ran = ran_johnson_new(seed, n, xi, lambda, gamma, delta);

  ran->type = RAN_JOHNSONSU;

  return ran;
}

/**
 * Draw from an unbounded Johnson variate.
 *
 * @param random is a generic random number generator.
 * @param xi is the localisation parameter.
 * @param lamdda is the scale parameter.
 * @param gamma and delta are the form parameters.
 *
 * @return the generated draw.
 */
double ran_johnsonSU(const Random *random,
		     double xi, double lambda, double gamma, double delta)
{  
   double z;

   z = (ran_normal_icdf(random,0,1) - gamma)/delta;
   return xi + lambda * (exp(z) - exp(-z))/2;
}

/**
 * Retrieve the parameters of an unbounded Johnson variate.
 *
 * @param random is an initialized random numbers generator.
 * @param xi contains on output the localisation parameter.
 * @param lamdda contains on ouput the scale parameter.
 * @param gamma and delta contain on output the form parameters.
 */
void ran_johnsonSU_get_parameters(const Random *random,
				  double *xi, double *lambda,
				  double *gamma, double *delta)
{
  ran_johnson_get_parameters(random, xi, lambda, gamma, delta);
}

/**
 * Set the parameters of an unbounded bounded Johnson random variate,
 *
 * @param random is an initialized random numbers generator.
 * @param xi is the localisation parameter.
 * @param lamdda is the scale parameter.
 * @param gamma and delta are the form parameters.
 */
void ran_johnsonSU_set_parameters(Random *random,
                                  double xi, double lambda, 
                                  double gamma, double delta)
{  
  ran_johnson_set_parameters(random, xi, lambda, gamma, delta);
}

/**
 * Draw from an unbounded Johnson variate unsing a specialized generator.
 *
 * @param random is a random number generator created specifically for
 *   some unbounded Johnson variate.
 * @param xi is the localisation parameter.
 * @param lamdda is the scale parameter.
 * @param gamma and delta are the form parameters.
 *
 * @return the generated draw.
 */
double ran_johnsonSU_get_val(const Random *random)
{
  double xi, lambda, gamma, delta; /* parameters of the Johnson
				      variate */

  ran_johnsonSU_get_parameters(random, &xi, &lambda, &gamma, &delta);

  return ran_johnsonSU(random, xi, lambda, gamma, delta);   
}

/**
 * Johnson SB: bounded form
 */

/**
 * Create a new bounded Johnson random variable generator with some
 * initial state.
 *
 * @param seed is the seed used to initialize the random numbers generator.
 * @param xi is the localisation parameter.
 * @param lamdda is the scale parameter.
 * @param gamma and delta are the form parameters.
 *
 * @return the newly allocated generator object.
 */
Random *ran_johnsonSB_new(unsigned long *seed, int n,
			  double xi, double lambda,
			  double gamma, double delta)
{
  Random *ran;

  ran = ran_johnson_new(seed, n, xi, lambda, gamma, delta);

  ran->type = RAN_JOHNSONSB;

  return ran;
}

/**
 * Draw from an bounded Johnson variate.
 *
 * @param random is a generic random number generator.
 * @param xi is the localisation parameter.
 * @param lamdda is the scale parameter.
 * @param gamma and delta are the form parameters.
 *
 * @return the generated draw.
 */
double ran_johnsonSB(const Random *random,
		     double xi, double lambda, double gamma, double delta)
{  
   double z;

   z = (ran_normal_icdf(random,0,1) - gamma)/delta;
   return xi + lambda * 1/(1 + exp(-z));
}

/**
 * Retrieve the parameters of a bounded Johnson variate.
 *
 * @param random is an initialized random numbers generator.
 * @param xi contains on output the localisation parameter.
 * @param lamdda contains on ouput the scale parameter.
 * @param gamma and delta contain on output the form parameters.
 */
void ran_johnsonSB_get_parameters(const Random *random,
				  double *xi, double *lambda,
				  double *gamma, double *delta)
{
  ran_johnson_get_parameters(random, xi, lambda, gamma, delta);
}

/**
 * Set the parameters of a bounded Johnson random variate,
 *
 * @param random is an initialized random numbers generator.
 * @param xi is the localisation parameter.
 * @param lamdda is the scale parameter.
 * @param gamma and delta are the form parameters.
 */
void ran_johnsonSB_set_parameters(Random *random,
				  double xi, double lambda, 
				  double gamma, double delta)
{
  ran_johnsonSB_set_parameters(random, xi, lambda, gamma, delta);
}

/**
 * Draw from an bounded Johnson variate unsing a specialized generator.
 *
 * @param random is a random number generator created specifically for
 *   some bounded Johnson variate.
 * @param xi is the localisation parameter.
 * @param lamdda is the scale parameter.
 * @param gamma and delta are the form parameters.
 *
 * @return the generated draw.
 */
double ran_johnsonSB_get_val(const Random *random)
{
  double xi, lambda, gamma, delta; /* parameters of the Johnson
				      variate */

  ran_johnsonSU_get_parameters(random, &xi, &lambda, &gamma, &delta);

  return ran_johnsonSB(random, xi, lambda, gamma, delta);   
}
