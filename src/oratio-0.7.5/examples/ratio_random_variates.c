#include <stdio.h>

#include <oratio/stat.h>
#include <oratio/ran.h>

/**
 * Simulate ratio of probabilities
 *
 * @param n is the number of simulations to perform.
 * @param r is a U(0,1) random numbers generator.
 */
void simulate_ratio(int n, Random *r)
{
  ST_TallyStore *num = st_tallystore_new();
  ST_TallyStore *den = st_tallystore_new();
  ST_TallyStore *ratio = st_tallystore_new();

  ST_TallyStore *boot_num = st_tallystore_new();
  ST_TallyStore *boot_den = st_tallystore_new();
  ST_Tally *boot_oups = st_tally_new();

  ST_TallyStore *boot_ratio_ok = st_tallystore_new();
  ST_TallyStore *boot_ratio_t = st_tallystore_new();
  ST_TallyStore *boot_ratio_oups = st_tallystore_new();

  double x, y;  /* random draws */
  double meanx, meany, ratio_mean, ratio_var;
  double boot_meanx, boot_meany, boot_ratio_mean, boot_ratio_var;
  double cov, boot_cov;
  double q;

  int i, k, m;
  int l, u;

  st_tallystore_set_name(num, "X");
  st_tallystore_set_name(den, "Y");
  st_tallystore_set_name(ratio, "ratio");

  st_tallystore_set_name(boot_ratio_ok, "Bootstrap Ratio OK");
  st_tallystore_set_name(boot_ratio_ok, "Bootstrap-t Ratio");
  st_tallystore_set_name(boot_ratio_oups, "Oups!");

  for (i = 0; i < n; i++) {
    /* Generate X following a normal of mean 2 and standard deviation 1 */
    x = ran_normal(r, 2, 1);
    /* Generate Y following a normal of mean 1 and standard deviation 2 */
    y = ran_normal(r, 1, 2);

    st_tallystore_add(num, x);
    st_tallystore_add(den, y);
    st_tallystore_add(ratio, x/y);
  }
  st_tallystore_report_basic(stdout, num);
  st_tallystore_report_ci_student(stdout, num, 0.95);
  st_tallystore_report_basic(stdout, den);
  st_tallystore_report_ci_student(stdout, den, 0.95);
  st_tallystore_report_basic(stdout, ratio);
  st_tallystore_report_ci_student(stdout, ratio, 0.95);

  meanx = st_tallystore_average(num);
  meany = st_tallystore_average(den);
  /* Compute the covariance between X and Y. We expect zero as the variables
   * has been indepedently generated. */
  cov = st_covariance(n, num->obs, 1, meanx, den->obs, 1, meany, 0, n);
  fprintf(stdout, "Covariance: %g\n", cov);
  /* Estimate the variance from the Delta theorem */
  ratio_mean = meanx/meany;
  ratio_var = (st_tallystore_variance(num)+
	       st_tallystore_variance(den)*ratio_mean*ratio_mean
	       -2*cov*ratio_mean)/(meany*meany);
  printf("Ratio average: %g\t Standard deviation: %g\n",
	 ratio_mean, sqrt(ratio_var));
  /* Delta confidence interval */
  q = st_gaussian_icdf(0.975);
  fprintf(stdout, "Delta confidence internal at level 0.95: [%g, %g]\n",
	  ratio_mean-q*sqrt(ratio_var/n), ratio_mean+q*sqrt(ratio_var/n));

  /* Recompute everything using bootstrap */
  m = 1000;
  for (i = 0; i < m; i++) {
    st_tallystore_init(boot_num);
    st_tallystore_init(boot_den);
    st_tally_init(boot_oups);

    for (k = 0; k < n; k++) {
      st_tallystore_add(boot_num, num->obs[ran_intrandom(r, n)]);
      st_tallystore_add(boot_den, den->obs[ran_intrandom(r, n)]);
      st_tally_add(boot_oups, ratio->obs[ran_intrandom(r, n)]);
    }

    st_tallystore_add(boot_ratio_ok,
		      st_tallystore_average(boot_num)/
		      st_tallystore_average(boot_den));

    boot_meanx = st_tallystore_average(boot_num);
    boot_meany = st_tallystore_average(boot_den);
    boot_ratio_mean = boot_meanx/boot_meany;
    boot_cov = st_covariance(n, boot_num->obs, 1, boot_meanx,
			     boot_den->obs, 1, boot_meany, 0, n);
    boot_ratio_var = (st_tallystore_variance(boot_num)+
		      st_tallystore_variance(boot_den)*
		      boot_ratio_mean*boot_ratio_mean
		      -2*boot_cov*boot_ratio_mean)/(boot_meany*boot_meany);
    st_tallystore_add(boot_ratio_t,
		      (boot_ratio_mean-ratio_mean)/sqrt(boot_ratio_var));

    st_tallystore_add(boot_ratio_oups, st_tally_average(boot_oups));
  }

  st_tallystore_report_basic(stdout, boot_ratio_oups);
  st_tallystore_report_basic(stdout, boot_ratio_ok);
  st_tallystore_report_ci_student(stdout, boot_ratio_ok, 0.95);

  l = ceil(m*(1-0.025));
  u = ceil(m*0.025);
  ds_sort_dheapsort(m, boot_ratio_ok->obs);
  ds_sort_dheapsort(m, boot_ratio_t->obs);
  printf("Bootstrap confidence interval: [%g, %g]\n",
	 2*ratio_mean-boot_ratio_ok->obs[l],
	 2*ratio_mean-boot_ratio_ok->obs[u]);
  printf("Bootstrap-t confidence interval: [%g, %g]\n",
	 ratio_mean-boot_ratio_t->obs[l]*sqrt(ratio_var),
	 ratio_mean-boot_ratio_t->obs[u]*sqrt(ratio_var));

  /* Bias estimation */
  printf("Estimated bias: %g\n",
	 st_tallystore_average(boot_ratio_ok)-ratio_mean);

  st_tallystore_free(num);
  st_tallystore_free(den);
  st_tallystore_free(ratio);

  st_tallystore_free(boot_num);
  st_tallystore_free(boot_den);
  st_tally_free(boot_oups);

  st_tallystore_free(boot_ratio_ok);
  st_tallystore_free(boot_ratio_t);
  st_tallystore_free(boot_ratio_oups);
}

int main()
{
  Random *r = ran_random();

  simulate_ratio(1000, r);
  printf("\n");
  simulate_ratio(10000, r);
  printf("\n");
  simulate_ratio(100000, r);
  printf("\n");
  simulate_ratio(1000000, r);
  printf("\n");

  ran_random_free(r);

  return 0;
}
