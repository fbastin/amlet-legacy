/**
 * Test with Oratio 0.7.1.
 * Compilation command: gcc normal_generation.c -lstat -lran -lctools
 */

#include <stdio.h>
#include <oratio/ctools.h>
#include <oratio/ran.h>

int main()
{
  Random *r;
  double x;

  long sec, ms;

  int i, n = 100000000;

  r = ran_normal_new(NULL, 0, 0.0, 1.0);

  ctools_timer_init();

  for (i = 0; i < n; i++)
    x = ran_gaussian_ar_laplace(r);
  ctools_timer(NULL, &sec, &ms);

  fprintf(stdout, "Time with acceptance-rejection method, for %d numbers: %ld.%ld seconds\n", n, sec, ms);

  ctools_timer_init();

  for (i = 0; i < n; i++)
    x = ran_gaussian_box_muller(r);
  ctools_timer(NULL, &sec, &ms);

  fprintf(stdout, "Time with Box-Muller method, for %d numbers: %ld.%ld seconds\n", n, sec, ms);

  ctools_timer_init();

  for (i = 0; i < n; i++)
    x = ran_gaussian_icdf(r);
  ctools_timer(NULL, &sec, &ms);

  fprintf(stdout, "Time with inversion method, for %d numbers: %ld.%ld seconds\n", n, sec, ms);

  ran_random_free(r);

  return 0;
}
