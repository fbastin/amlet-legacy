/**
 * io.h
 *
 * OPHELIA Library
 * Copyright (C) 1999-2010 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

/**
 * Output Functions
 */

#ifndef __NLP2_IO_H__
#define __NLP2_IO_H__

#include "problem.h"
#include "outer.h"
#include "btr_struct.h"
#include "bfgs_struct.h"

#ifdef __cplusplus
extern "C" {
#endif

  enum {
    SILENT = 0,
    MODERATE_VERBOSE,
    HIGHLY_VERBOSE
  };

  /**
   * Print the problem name in the output file.
   *
   * @param output is the file where to print name.
   * @param p is the problem under consideration.
   */
  void op_problem_print_name(FILE *output, OP_Problem *p);

  /**
   * Outer iterations
   */

  void nlp_outer_print_iteration(FILE *output, nlp_outer *no, void *data);

  /**
   * BFGS I/O functions
   */

  void bfgs_print_iteration(FILE *output, BFGS *bf, void *data);

  void bfgs_print_results(FILE *output, BFGS *bf, void *data);

  /**
   * BTR I/O functions
   */

  void btr_print_results(FILE *output, BTR *btr, void *data);

  void btr_print_iteration(FILE *output, BTR *btr, void *data);

  void btr_print_gnuplot(FILE *output, BTR *btr, void *data);

  void btr_print_iteration_latex(FILE *output, BTR *btr, void *data);

#ifdef __cplusplus
}
#endif

#endif /* __NLP2_IO_H__ */
