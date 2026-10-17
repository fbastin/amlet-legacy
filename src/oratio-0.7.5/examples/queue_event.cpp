/**
 * M/M/1 queue simulation.
 *
 * Tested with Oratio 0.7.2.
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
  void simulates(int n);

  void empty();
};

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

void customer_free(void *c)
{
  delete((Customer *)c);
}

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

Queue::Queue(double lambda, double mu)
{
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

Queue::~Queue()
{
  ran_random_free(gen_arr);
  ran_random_free(gen_serv);
  /* to improve, but I have to be careful with the mix of C and C++ */
  ds_list_free(wait_list, customer_free);
  ds_list_free(serv_list, customer_free);
  delete(sim);
  st_tally_free(cust_waits);
  delete (tot_wait);
}

void Queue::empty()
{
  ds_list_free(wait_list, customer_free);
  ds_list_free(serv_list, customer_free);
  wait_list = ds_list_new();
  serv_list = ds_list_new();
  st_tally_init(cust_waits);
  tot_wait->init();
}

void Queue::simulate(double time_horizon)
{
  sim->init();
  empty();

  EndOfSim *eos = new EndOfSim(sim);
  eos->schedule(time_horizon);

  Arrival *a = new Arrival(sim, this);
  a->schedule(ran_random_next(gen_arr));

  sim->start();
}

void Queue::simulates(int n)
{
  for (int i = 0; i < n; i++) {
    simulate(1000.0);


    cout << "Simulation " << i << endl;
    st_tally_report_basic(stdout, cust_waits);
    st_tally_report_ci_student(stdout, cust_waits, 0.95);
    tot_wait->report(stdout);
  }
}

int main() {
  Queue *q = new Queue(1.2, 2.0);

  q->simulates(10);

  delete(q);

  return 0;
}
