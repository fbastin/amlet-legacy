/**
 * Test with Oratio 0.7.1.
 * Compilation command: gcc normal_ar.c -lstat -lran -lctools
 */

#include <stdio.h>
#include <oratio/ctools.h>
#include <oratio/ran.h>
#include <oratio/stat.h>

int main()
{
  Random *r;
  double x;

  ST_Tally *st = st_tally_new();

  long sec, ms;

  int i, n = 1000000;

  r = ran_normal_new(NULL, 0, 0.0, 1.0);

  ctools_timer_init();

  for (i = 0; i < n; i++) {
    x = ran_gaussian_ar_laplace(r);
    st_tally_add(st, x);
  }

  fprintf(stdout, "Mean: %g; variance: %g\n", st_tally_average(st),
	  st_tally_variance(st));

  st_tally_free(st);

  ran_random_free(r);

  return 0;
}
