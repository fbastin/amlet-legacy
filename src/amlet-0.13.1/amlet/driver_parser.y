/**
 * driver_parser.y
 *
 * AMLET (Another Mixed Logit Estimation Tool)
 * Copyright (C) 2002-2010 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

%{

#include <stdio.h>
#include <string.h>

#include "amlet_driver.h"

#include <oratio/ctools.h>

  extern FILE *yyin;
  extern void yyrestart(FILE *input_file);
  int  yyerror(char *s);

  int lineno;
  AMLET_Driver *dr;

  driver_var *stack;
  int max_stack;
  int n = 0;

  int type;  /* not very a babe this variable */

  void push_stack(driver_var drv)
    {
      if (n >= max_stack) {
	max_stack *= 2;
	stack = _reallocs(stack, max_stack, driver_var);
      }

      stack[n++] = drv;
    }

  %}

%union {
  double  real;      /* real value */
  int     integer;   /* integer value */
  char*   string;    /* string buffer */
  void*   ptr;       /* just for reading facility */
}

%token <integer>  INTEGER IDENTIFIER
%token <real>     REAL
%token <string>   VARIABLE
%token <integer>  POS
%token <string>   STRING

%token <integer>  DISTRIBUTION
%token <integer>  OPTMETHOD
%token <integer>  HESSIANUPDATE
%token <integer>  SAMPLINGMETHOD
%token <integer>  RANDOMIZATIONMETHOD
%token <integer>  RANGENERATOR
%token <integer>  FORMAT_DATA
%token <integer>  NONPARAMETRIC

%type <integer>  iexp;  /* integer expression */
%type <real>     dexp;  /* real expression */
%type <ptr>  identifier;

%type <ptr>      var;

%token LBRACE   /* '{' */
%token RBRACE   /* '}' */

%token IDENTIFIER  /* identifier of a reserved word */

/* commands */
%token PRINT
%token LOAD
%token GENERATE
%token INPUT
%token OUTPUT

%%

source:
        instructions {driver_check(dr);}
        ;

instructions:
	instructions instruction {;}
      |
      ;

instruction: inst ';'
      | ';'   /* we allow empty instructions */
      ;

inst: affectation
      | command
      ;

command:
         generate   {dr->generation = 1;}
         | load
	 | input
	 | output
         | PRINT
         ;

generate:
	 GENERATE
           {
	     driver_set_generation_process(dr, NULL, dr->var[NALT].i);
	   }
         | GENERATE STRING
           {
	     driver_set_generation_process(dr, $2, dr->var[NALT].i);
	   }
         | GENERATE STRING STRING
           {
	     driver_set_generation_process(dr, $2, dr->var[NALT].i);
	     fprintf(stderr, "** WARNING ** this kind of panel generation is to be implemented\n");
           }
         ;

load:
      LOAD STRING
        {
	  if (dr->generation)
	    c_warning("observations file specified after a generating process ghas been given.\n");
	  dr->generation = 0;
	  dr->obsname = $2;
	}
      ;

input: INPUT STRING
         {
	   FILE *f2;

	   f2 = fopen($2, "r");
         }
       ;

output:
      OUTPUT STRING
        {
          dr->output = $2;
	}
      ;

affectation:
             var '=' iexp {;}
             | var '=' dexp {;}
             | var '=' STRING {free($3);}
             | var '=' set
               {
		 int i;

		 if (type == INTEGER) {
		   int *array;
		   array = _mallocs(n, int);
		   for (i = 0; i < n; i++) array[i] = stack[i].d.i;
		   (*((int **)$1)) = array;
		 }
		 else if (type == REAL) {
		   double *array;
		   array = _mallocs(n, double);
		   for (i = 0; i < n; i++) {
		     if (stack[i].type == INTEGER)
		       array[i] = stack[i].d.i;
		     else
		       array[i] = stack[i].d.d;
		   }
		   (*((double **)$1)) = array;
		 }
		 else
		   $1 = NULL;
		 n = 0;
	       }
             | identifier '=' iexp
               {
		 (*((driver_data *)$1)).i = $3;
               }
             | identifier '=' dexp {(*((driver_data *)$1)).d = $3;}
             | identifier '=' STRING {(*((driver_data *)$1)).s = $3;}
             | identifier '=' set
               {
		 int i;
		 driver_array *a = (driver_array *)$1;

		 if (a->array) free(a->array);

		 a->n = n;
		 a->array = (!n) ? NULL: _mallocs(n, driver_var);
		 for (i = 0; i < n; i++) {
		   a->array[i].type = stack[i].type;
		   a->array[i].d = stack[i].d;
		 }

		 /* We return to the beginning of the stack */
		 n = 0;
               }
             ;

identifier: IDENTIFIER
            {
	      int pos = reserved_words[$1].token;

	      if ( yylval.integer <
		   sizeof(reserved_words)/sizeof(reserved_words[0])
		   - END_ARRAY_WORDS - END_OTHERS )
		/* Simple variable */ {
		$$ = (dr->var + pos);
	      }
	      else if ( yylval.integer <
			sizeof(reserved_words)/sizeof(reserved_words[0])
			- END_OTHERS)
		/* Array variable */
		$$ = (dr->array + pos);
	      else
		driver_get_address($1);
	    }
            | IDENTIFIER POS { $$ = NULL; }
            ;

var: VARIABLE {;}
     | VARIABLE POS
       {
	 int r, pos;

	 /* Determine */
	 r = c_symbol_search($1, arrays,
			     sizeof(arrays)/sizeof(arrays[0]), &pos);

	 if (r < 0 || r >= END_DBL_ARRAY_WORDS) {
	   c_error("unrecognized array element at line %d\n", lineno);
	 }

	 if (!driver_check_array_pos(dr, r, $2)) {
	   c_error("invalid array position in line %d\n", lineno);
	 }

	 $$ = driver_array_address(dr, r, $2, &type);
       }
     ;

iexp:
      INTEGER
      | DISTRIBUTION
      | NONPARAMETRIC
      | FORMAT_DATA
      | HESSIANUPDATE
      | OPTMETHOD
      | RANGENERATOR
      | SAMPLINGMETHOD
      | RANDOMIZATIONMETHOD
      ;

dexp:
      REAL
      ;

set: LBRACE values RBRACE {;}
     ;

string:
        STRING
        {
          driver_var drv;
	  drv.type = STRING;
	  drv.d.s = $1;
	  push_stack(drv);
        }
        ;

values:
        value
        | values ',' value
        | values ',' string
        |
        ;

value:
	iexp
          {
	    driver_var drv;
	    drv.d.i = $1;
	    drv.type = INTEGER;
	    push_stack(drv);
          }
        | dexp
          {
            driver_var drv;
	    drv.d.d = $1;
	    drv.type = REAL;
	    push_stack(drv);
          }
        ;

%%

void amlet_parser_init(AMLET_Driver *dri)
{
  lineno = 1;

  dr = dri;

  stack = _mallocs(32, driver_var);
  max_stack = 32;
}

void amlet_parser_clean()
{
  /* quid of the string elemets? */
  free(stack);
}

int yyerror(char *s) {
  fprintf(stderr, "%s at line %d\n", s, lineno);
  return 0;
}
