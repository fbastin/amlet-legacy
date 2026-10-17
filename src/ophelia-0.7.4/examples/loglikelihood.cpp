/**
 * Compilation: g++ loglikelihood.cpp -lntools -lran -lstat -lophelia
 */

#include <iostream>
#include <fstream>

using namespace std;

#include <oratio/oratio.h>

#include <ophelia/btr_struct.h>
#include <ophelia/ophelia.hpp>

typedef struct _queue {
  ST_TallyStore *arrivals;
  ST_Obs *obs;
} queue;

inline void obs_component(ST_Obs *obs, double x, double lambda)
{
  double e = 1.0e-8;
  obs->prob = st_exponential_pdf(x, lambda);
  obs->grad[0] = st_exponential_pdf_llderivative(x, lambda);
}

void observations(ST_TallyStore *arrivals, ST_Obs *obs, double lambda)
{
  for (int i = 0; i < st_tallystore_nobs(arrivals); i++)
    obs_component(obs+i, arrivals->obs[i], lambda);
}

void gradient(double *x, int n, void *data, double *grad)
{
  queue *q = (queue *)data;
  ST_TallyStore *arrivals = q->arrivals;
  ST_Obs *obs = q->obs;

  st_loglikelihood_grad(st_tallystore_nobs(arrivals), 1, obs, grad, 1);

  *grad = -(*grad);
}

double objective(double *x, int n, void *data)
{
  queue *q = (queue *)data;
  ST_TallyStore *arrivals = q->arrivals;
  ST_Obs *obs = q->obs;

  int nobs = st_tallystore_nobs(arrivals);

  observations(arrivals, obs, *x);

  return -st_loglikelihood(nobs, obs);
}

double fg(double *x, int n, void *data, double *grad)
{
  gradient(x, n, data, grad);
  return objective(x, n, data);
}
/**
 * Generate a Poisson arrival process.
 *
 * @lambda is the parameter of the arrival process; its inverse gives the
 *   mean time between two arrivals.
 */
void generate_poisson(double lambda, ST_TallyStore *arrivals)
{
  ofstream out;
  Random *r = ran_exponential_new(NULL, 0, 1.0);
  double time;

  out.open("poisson_arrival.dat");
  if (out.is_open()) {
    for (int i = 0; i < 100000; i++) {
      time = ran_random_next(r);
      st_tallystore_add(arrivals, time);
      out << time << endl;
    }
    out.close();
  }

  ran_random_free(r);
}

int main()
{
  ST_TallyStore *arrivals = st_tallystore_new();
  st_tallystore_set_name(arrivals, "Arrivals");
  generate_poisson(1.0, arrivals);
  st_tallystore_report_basic(stdout, arrivals);

  double lambda;
  double **H;

  cout << "Maximum likelihood estimator of the arrival rate: "
       << st_exponential_ml_estimator(st_tallystore_nobs(arrivals),
				      arrivals->obs, 1, &lambda)
       << endl;

  // Solve the maximum likelihood estimation problem

  // Initial lambda
  lambda = 0.5;

  queue q;
  q.arrivals = arrivals;
  q.obs = new ST_Obs[st_tallystore_nobs(arrivals)];
  for (int i = 0; i < st_tallystore_nobs(arrivals); i++)
    q.obs[i].grad = new double;

  BTR *b = op_btr_new();

  op_btr_init(b, 1, 0);

  /* Starting point */
  b->x[0] = lambda;

  b->printer = btr_print_iteration;

  H = nt_matrix_new(1,1);
  nlp_btr(b, (C_GENERIC)objective, (C_GENERIC)gradient,
  	  &q, NULL, H, NULL);
  /*
    nlp_btr(b, (C_GENERIC)objective, NULL,
    &q, NULL, H, NULL);
  */

  // Fill the observations

  lambda = b->x[0];
  cout << "Maximum likelihood estimator of the arrival rate: "
       << lambda << endl;

  op_btr_free(b);
  nt_matrix_free(H);

  for (int i = 0; i < st_tallystore_nobs(q.arrivals); i++) delete(q.obs[i].grad);
  delete(q.obs);
  st_tallystore_free(arrivals);

  return 0;
}
