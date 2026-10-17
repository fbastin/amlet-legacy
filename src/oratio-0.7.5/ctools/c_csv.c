/**
 * c_csv.c
 *
 * Copyright (C) 2000-2008 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include <stdlib.h>
#include <string.h>

#include <oratio/ctools.h>

static char comments='#';

#define C_CSV_DEFAULT_SEPARATOR ';'

void *c_csv_malloc(size_t s)
{
  void *p;

  if (!(p = malloc(s))) {
    fprintf(stderr,"Memory leaks for CSV library\nI will take some rest until you give me a more pleasant workspace.\nBye!\n");
    exit(1);
  }

  return p;
}

/**
 * @param separator is the separator that will be used between
 *   entries. Separator must be a pointer to a single character, or set
 *   to NULL, in which C_CSV_DEFAULT_SEPARATOR will be used.
 */
C_CSV *c_csv_new(char separator)
{
  C_CSV *csv;

  csv = (C_CSV *)c_csv_malloc(sizeof(C_CSV));
  csv->blanks = strdup(" \t");
  csv->separator = (separator != '\0') ? separator : C_CSV_DEFAULT_SEPARATOR;
  csv->comments = strdup("#");

  return csv;
}

/**
 * Free the memory associated to a C_CSV object.
 */
void c_csv_free(C_CSV *csv)
{
  if (csv) {
    free(csv->blanks);
    free(csv->comments);
    free(csv);
  }
}

#define isComment(csv, line) (line == strstr(line, csv->comments))

int is_comment(char *line) {
  if (line[0] == comments)
    return 1;
  else
    return 0;
}

int count_nb_items(C_CSV *csv, char *line)
{
  int i = 0;
  int nb_items = 1;
  int quote = 0;

  while (line[i] != '\n' && line[i] != EOF) {
    if (line[i] == '\"') {
      if (quote) quote = 0; else quote = 1;
    }
    else if (line[i] == csv->separator) {
      if (!quote) nb_items++;
    }
    i++;
  }

  if (nb_items == 1) {
    i = 0;
    while (strchr(" \t", line[i])) i++;
    if (line[i] == '\n') nb_items = 0;
  }

  return nb_items;
}

char *read_csv_item(C_CSV *csv, char *line, int *pos)
{
  int eoln = 0;
  int quote = 0;
  int first = *pos, pos2 = *pos;
  char *buff;

  if (csv->separator != '\t' && csv->separator != ' ')
    /* Skip spaces and tabulations */
    while (strchr(" \t", line[*pos])) (*pos)++;

  while ((!eoln && line[*pos] != csv->separator) ||
	 (line[*pos] == csv->separator && quote)) {
    switch (line[*pos]) {
    case '\"':
      if (line[*pos+1] == '\"')
	line[pos2++] = '\"';
      else {
	if (quote)
	  quote = 0;
	else
	  quote = 1;
      }
      break;
    case EOF:
    case '\n':
      eoln = 1;
      break;
    case '\r':
      if (line[*pos+1] == '\n') {
	eoln = 1;
	break;
      }
    default:
      line[pos2++] = line[*pos];
    }
    (*pos)++;
  }

  buff = (char *)malloc((pos2-first+1)*sizeof(char));
  memcpy(buff, &(line[first]), pos2-first);
  buff[pos2-first] = '\0';

  return buff;
}

/**
 * @param length is a pointer to an integer whose value is on output
 *   the number of items contained in the current line.
 *
 * @return an array of characters strings containing all the items of
 *   the line; this array is newly allocated, and consequently, should
 *   be correctly deallocated after use.
 */
char **c_csv_read_line(C_CSV *csv, FILE *input, int *length)
{
  char *line, **items;
  int i, pos = 0;

  items = NULL;

  /* Read the line */
  if ((line = c_file_get_line(input))) {
    if (line) {
      if (line[0] == '\n' || is_comment(line) || line[0] == '\0') {
	free(line);
	line = NULL;
      }
    }
  }
  *length = 0;

  if (line) {
    /* Treat the line */
    *length = count_nb_items(csv, line);
    if (*length) {
      items = (char **)malloc(1+*length*sizeof(char *));

      for (i = 0; i < *length; i++) {
	items[i] = read_csv_item(csv, line, &pos);
	pos++; /* on passe le séparateur */
      }
      items[*length] = NULL;
    }
    free(line);
  }

  return items;
}

void write_csv_item(FILE *output, char *item, char separator)
{
  int quote = 0, j=0;

  if (strchr(item, separator) ||
      strchr(" \t", item[0]) ||
      strchr(" \t", item[strlen(item)-1])) quote = 1;

  if (quote) fputc('\"', output);

  while (item[j] != '\0') {
    if (item[j] == '\"')
      fputs ("\"\"", output);
    else
      fputc (item[j], output);
    j++;
  }

  if (quote) fputc('\"', output);
}

void c_csv_write_line(C_CSV *csv, FILE *output,
		      char **items, int nb_items)
{
  int i;

  for (i = 0; i < nb_items-1; i++) {
    if (items[i]) write_csv_item(output, items[i], csv->separator);
    fputc(csv->separator, output);
  }
  write_csv_item(output, items[nb_items-1], csv->separator);

  fprintf(output, "\r\n");
}
