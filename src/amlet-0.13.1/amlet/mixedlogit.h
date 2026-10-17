#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Compute the simulation bias.
   *
   * Ref.: F. Bastin, C. Cirillo, and Ph. L. Toint, An adaptive Monte Carlo
   * algorithm for computing mixed logit estimators, Computational Management
   * Science, 3(1):55--79 (2006).
   * Special thanks to Adrien Bonneu, INSA, Toulouse. First bias correction.
   *
   * @param popsize is the population size
   * @param accuracy is the estimated precision.
   * @param alpha is the significance level used for the accuracy estimation.
   *
   * @return the bias of simulation
   */
  double simulation_bias(int popsize, double alpha, double accuracy);

#ifdef __cplusplus
}
#endif
