#include <stdio.h>

#include <oratio/ran.h>
#include <oratio/stat.h>

int main()
{
  Random *r;
  ST_TallyStore *set = st_tallystore_new();
  ST_TallyStore *maxest = st_tallystore_new();

  double theta = 1.0;
  double x, max;
  int count = 0;

  int i, j;
  int n = 5000, m = 1000;

  r = ran_random();

  /* Generate observations */
  for (i = 0; i < n; i++) st_tallystore_add(set, ran_random_next(r)*theta);

  for (i = 0; i < m; i++) {
    max = 0.0;
    for (j = 0; j < n; j++) {
      x = set->obs[ran_intrandom(r, n)];
      if (x > max) max = x;
    }
    st_tallystore_add(maxest, max);
    if (max == set->tally->max_value) count++;
  }

  st_tallystore_report_basic(stdout, set);
  st_tallystore_report_basic(stdout, maxest);
  printf("Frequence of maximal observation report: %g\n",((double)count)/m);

  ran_random_free(r);
  st_tallystore_free(set);
  st_tallystore_free(maxest);

  return 0;
}
