/**
 * parabola_1.c
 *
 * OPHELIA Library
 * Copyright (C) 2004-2015 Fabian Bastin <fabian.bastin@cirrelt.ca>
 *
 * Licensed under the MIT License
 */

/**
 * Evaluate the function x²+4x-10.
 *
 * @param x is the point where the function has to be evaluated.
 * @param data is unused.
 *
 * @return the value of x²+4x-10.
 */
double op_parabola_1(double x, void *data)
{
  return x*(x+4)-10;
}
