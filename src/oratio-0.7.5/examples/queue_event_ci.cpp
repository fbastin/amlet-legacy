/**
 * M/M/1 queue simulation
 * This example provides confidence intervals on the produced results.
 *
 * Tested with oratio 0.7.2
 * Compilation: g++ queue_event.cpp -lorsim -lran -lstat -lntools -ldstruct
 */

#include <iostream>
using namespace std;

#include <stdio.h>

#include <math.h>

#include <oratio/dstruct.h>
#include <oratio/ran.h>
#include <oratio/sim_accumulate.hpp>
#include <oratio/sim.hpp>
#include <oratio/st_tally.h>

/* System to simulate */
class Queue
{
  ESim *sim;

  double lambda, mu;
  Random *gen_arr;
  Random *gen_serv;

  DS_List *wait_list;
  DS_List *serv_list;

public:

  ST_Tally *cust_waits;
  ESim_Accumulate *tot_wait;

  class Customer {
    double arrivTime, servTime;
    int id;

  public:
    Customer(double ar, double se) {
      arrivTime = ar;
      servTime = se;
    };

    friend class Queue;
  };

  class Arrival : public ESim_Event {
    Queue *q;

  public:

    Arrival(ESim* sim, Queue *queue) : ESim_Event(sim) {q = queue;};

    void actions();
  };

  class Departure : public ESim_Event {
    Queue *q;

  public:

    Departure(ESim* sim, Queue *queue) : ESim_Event(sim) {q = queue;};

    void actions();
  };

  class EndOfSim : public ESim_Event {    

  public:
    EndOfSim(ESim* sim) : ESim_Event(sim) {};

    void actions() {
      sim->stop();
    };
  };

public:

  Queue(double lambda, double mu);
  virtual ~Queue();

  void simulate(double time_horizon);
};

void Queue::Arrival::actions()
{
  Arrival *a = new Arrival(sim, q);

  a->schedule(ran_random_next(q->gen_arr));

  Customer *c = new Customer(sim->clock, ran_random_next(q->gen_serv));

  if (ds_list_length(q->serv_list)) {
    // Must join the queue.
    q->wait_list = ds_list_append(q->wait_list, c);
    q->tot_wait->update(ds_list_length(q->wait_list));
  }
  else {
    // Starts service.
    st_tally_add(q->cust_waits, 0.0);
    q->serv_list = ds_list_append(q->serv_list, c);
    Departure *d = new Departure(sim, q);
    d->schedule(c->servTime);
  }
}

void Queue::Departure::actions()
{
  Customer *c = (Customer *)ds_list_get_data(ds_list_first(q->serv_list));
  q->serv_list = ds_list_pop(q->serv_list);
  delete (c);

  if (ds_list_length(q->wait_list)) {
    q->wait_list = ds_list_first(q->wait_list);
    c = (Customer *)ds_list_get_data(q->wait_list);
    q->wait_list = ds_list_pop(q->wait_list);
    
    q->tot_wait->update(ds_list_length(q->wait_list));
    st_tally_add(q->cust_waits, q->sim->clock - c->arrivTime);
    q->serv_list = ds_list_append(q->serv_list, c);

    Departure *d = new Departure(sim, q);
    d->schedule(c->servTime);
  }
}

Queue::Queue(double lambda, double mu) {
  sim = new (ESim);

  this->lambda = lambda;
  this->mu = mu;

  gen_arr = ran_exponential_new(NULL, 0, lambda);
  gen_serv = ran_exponential_new(NULL, 0, mu);

  wait_list = ds_list_new();
  serv_list = ds_list_new();

  cust_waits = st_tally_new();
  tot_wait = new ESim_Accumulate(sim);
}

Queue::~Queue() {
  ran_random_free(gen_arr);
  ran_random_free(gen_serv);
  /* to improve, but I have to be careful with the mix of C and C++ */
  ds_list_free(wait_list, NULL);
  ds_list_free(serv_list, NULL);
  delete(sim);
  st_tally_free(cust_waits);
  delete (tot_wait);
}

void Queue::simulate(double time_horizon) {
  sim->init();

  EndOfSim *eos = new EndOfSim(sim);
  eos->schedule(time_horizon);

  Arrival *a = new Arrival(sim, this);
  a->schedule(ran_random_next(gen_arr));

  sim->start();
}

int main() {
  Queue *q = new Queue(1.99, 2.0);

  q->simulate(200.0);

  cout << "Temps d'attente moyen: " << st_tally_average(q->cust_waits) << endl;
  cout << "Taille moyenne de la file d'attente: " << q->tot_wait->average() << endl;
  st_tally_report_basic(stdout, q->cust_waits);
  st_tally_report_ci_normal(stdout, q->cust_waits, 0.95);
  st_tally_report_ci_student(stdout, q->cust_waits, 0.95);
  q->tot_wait->report(stdout);

  delete(q);

  return 0;
}
