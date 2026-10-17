/**
 * io.c
 *
 * ORATIO Library
 * Copyright (C) 2002-2009 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#include <oratio/ntools.h>
#include <oratio/ctools.h>

#include <string.h>

/**
 * Matrix input and output
 */

/**
 * Read the matrix A in the file in.
 *
 * @param in is the input file, containing the matrix to read.
 * @param name contains on output the matrix name.
 * @param A is on input a pointer to a NULL matrix or a matrix that will be
 *   destroyed, and on output a pointer to a newly created containing the
 *   read elements.
 * @param m is the number of rows of the matrix.
 * @param n is the number of columns of the matrix.
 */
int nt_matrix_read(FILE *in, char **name,
		   double ***A, int *m, int *n)
{
  char *line, *buf, *token;
  char a, b;
  int i, j, size = 0, smax;
  int err = 0;
  double tmp;
  double *values, *elem;
  int finished = 0;

  line = c_file_get_line(in);
  buf = _mallocs(strlen(line), char);

  /* Read the matrix declaration. */
  size = sscanf(line, "%s %c %c", buf, &a, &b);
  if (size != 3) {
    free(buf);
    free(line);
    return 1;
  }
  if (a != '=' || b != '[') {
    free(buf);
    free(line);
    return 1;
  }
  if (*name)
    *name = c_realloc(*name, strlen(buf)*sizeof(char));
  else
    *name = _mallocs(strlen(buf), char);
  strcpy(*name, buf);

  /* To improve */
  if (*A) nt_matrix_free(*A);
  *A = NULL;

  /* Read the matrix. */
  *m = *n = 0;
  smax = 1024;
  values = _mallocs(smax,double);

  /* The first line allows us to know the number of columns. */
  strcpy(buf, line);
  strtok(buf, " \t\n");
  for (i = 0; i < 3; i++) token = strtok(NULL, " \t\n");

  if (!sscanf(token, "%lg", &tmp)) {
    free(line);
    line = c_file_get_line(in);
    buf = c_realloc(buf, strlen(line)*sizeof(char));
  }
  size = 0;
  do {
    if (!sscanf(token, "%lg", &tmp)) {
      finished = 1;
      if (*token != ']') err = 1;
      break;
    }
    if (size >= smax) {
      smax <<= 1;
      values = c_realloc(values, smax*sizeof(double));
    }
    values[size++] = tmp;
    (*n)++;
  } while ((token = strtok(NULL, " \t\n")));
  if (size) *m = 1;

  /* Read the following lines. */
  while (!finished) {
    line = c_file_get_line(in);
    if (!line) {
      finished = err = 2;
      break;
    }
    strcpy(buf, line);
    token = strtok(buf, " \t\n");
    if (sscanf(line, "%lg", &tmp)) {
      (*m)++;
      /* Read the entries. */
      for (i = 0; i < *n; i++) {
	if (!sscanf(token, "%lg", &tmp)) {
	  finished = 1;
	  err = 3;
	}
	else {
	  if (size >= smax) {
	    smax <<= 1;
	    values = c_realloc(values, smax*sizeof(double));
	  }
	  values[size++] = tmp;
	}
	token = strtok(NULL, " \t\n");
      }
    }
    if (token) {
      finished = 1;
      if (*token != ']') err = 1;
    }
  }

  free(line);
  free(buf);

  if (!err) {
    /* Write the matrix in memory. */
    elem = values;
    *A = nt_matrix_new(*m, *n);
    for (i = 0; i < *m; i++)
      for (j = 0; j < *n; j++) {
	(*A)[i][j] = *elem;
	elem++;
      }
  }

  free(values);
  return err;
}

/**
 * Print a vector on a file.
 *
 * @param out is the file where the vector has to be written.
 * @param x is the vector to print.
 * @param n is dimension of x.
 */
void nt_vect_print(FILE *out, char *name, double *x, int n)
{
  int i;

  if (name) fprintf(out, "%s = ", name);
  fprintf(out, "[");
  for (i = 0; i < n; i++) fprintf(out, " %10.6g", x[i]);
  fprintf(out, " ]\n");
}

/**
 * Print the matrix A in the file out.
 *
 * @param out is the file where the matrix has to be written.
 * @param m is the number of rows of the matrix.
 * @param n is the number of columns of the matrix.
 */
void nt_matrix_print(FILE *out, char *name,
		     double **A, int m, int n)
{
  int i, j;
  int len;

  len = strlen(name);

  if (name) fprintf(out, "%s = [", name);
  for (i = 0; i < m; i++) {
    if (i) for (j = 0; j < len+4; j++) fprintf(out, " ");
    for (j = 0; j < n; j++) fprintf(out, " %12.8g", A[i][j]);
    if (i == m-1 && name)
      fprintf(out, " ]\n");
    else
      fprintf(out, "\n");
  }
}

/**
 *
 * @param out is the file where the matrix has to be written.
 * @param m is the number of rows of the matrix.
 * @param n is the number of columns of the matrix.
 */
void nt_matrix_print2(FILE *out, char *name,
		      double *A, int m, int n,
		      int lda)
{
  int i, j;
  int len;

  len = strlen(name);

#define a_ref(a1,a2) A[(a2)*(a_dim1)+a1]

  if (name) fprintf(out, "%s = [", name);
  for (i = 0; i < m; i++) {
    if (i) for (j = 0; j < len+1; j++) fprintf(out, " ");
    for (j = 0; j < n; j++)
      fprintf(out, " %12.8g", A[i*lda+j]);
    if (i == m-1 && name)
      fprintf(out, " ]\n");
    else
      fprintf(out, "\n");
  }
}

/**
 * Sets input/output.
 */

/**
 * Print a set into a file.
 *
 * @param out is the file where the set has to be written.
 * @param n is the set size.
 * @param s is the set to print.
 * @param lds is the increment of s.
 */
void nt_set_print(FILE *out, int n, double *s, int lds)
{
  int i;

  fprintf(out, "{ ");
  for (i = 0; i < n-1; i++) fprintf(out, "%10.6g, ", s[i*lds]);
  if (n) fprintf(out, "%10.6g ", s[(n-1)*lds]);
  fprintf(out, "}\n");
}

/**
 * Print a set of integers into a file.
 *
 * @param out is the file where the set has to be written.
 * @param n is the set size.
 * @param s is the set to print.
 * @param lds is the increment of s.
 */
void nt_set_iprint(FILE *out, int n, int *s, int lds)
{
  int i;

  fprintf(out, "{ ");
  for (i = 0; i < n-1; i++) fprintf(out, "%d, ", s[i*lds]);
  if (n) fprintf(out, "%d ", s[(n-1)*lds]);
  fprintf(out, "}\n");
}
