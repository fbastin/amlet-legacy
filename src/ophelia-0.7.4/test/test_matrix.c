/**
 * test_matrix.c
 *
 * OPHELIA Library
 * Copyright (C) 2011 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include <oratio/oratio.h>
#include <ophelia/matrix.h>
#include <ophelia/clapack.h>

void test_matrix()
{
  double **H;

  H = nt_matrix_new(3,3);

  H[0][0] = 1;
  H[0][1] = 2;
  H[0][2] = 3;
  H[1][0] = 2;
  H[1][1] = 5;
  H[1][2] = 6;
  H[2][0] = 3;
  H[2][1] = 6;
  H[2][2] = 2;

  nt_matrix_print(stdout, "H", H, 3, 3);
  op_matrix_syin_inverse(CblasColMajor, CblasUpper, *H, 3);
  op_matrix_copy_lower(CblasRowMajor, 3, *H, 3);
  nt_matrix_print(stdout, "H^-1", H, 3, 3);
  H[0][0] = 1;
  H[0][1] = 2;
  H[0][2] = 3;
  H[1][0] = 2;
  H[1][1] = 5;
  H[1][2] = 6;
  H[2][0] = 3;
  H[2][1] = 6;
  H[2][2] = 2;

  nt_matrix_print(stdout, "H", H, 3, 3);
  op_matrix_sy_inverse(CblasRowMajor, CblasUpper, *H, 3);
  nt_matrix_print(stdout, "H^-1", H, 3, 3);

  H[0][0] = 3;
  H[0][1] = 1;
  H[0][2] = 1;
  H[1][0] = 1;
  H[1][1] = 4;
  H[1][2] = 2;
  H[2][0] = 1;
  H[2][1] = 2;
  H[2][2] = 5;
  nt_matrix_print(stdout, "H", H, 3, 3);
  op_matrix_sy_inverse(CblasRowMajor, CblasUpper, *H, 3);
  nt_matrix_print(stdout, "H^-1", H, 3, 3);

  nt_matrix_free(H);
}
