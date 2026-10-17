/**
 * driver_scanner.l
 *
 * AMLET (Another Mixed Logit Estimation Tool)
 * Copyright (C) 2002-2010 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

/**
 * Scanner for the driver file
 */

/* The "incl" state is used for picking up the name of an include file.
 * The implementation follows the flex manpage.                         */
%x incl

%{
#include "driver_parser.h"

#include <string.h>

#include "amlet_data.h"
#include "amlet_driver.h"

#include <oratio/ctools.h>

#define MAX_INCLUDE_DEPTH 10
YY_BUFFER_STATE include_stack[MAX_INCLUDE_DEPTH];
int include_stack_ptr = 0;

  /* commands */
  static const c_symbol commands[] = {
    {"print", PRINT},
    {"load", LOAD},
    {"generate", GENERATE},
    {"output", OUTPUT}
  };

  /* functions declarations */
  int c_symbol_search(char *s, const c_symbol *symbols, int n, int *pos);
  int isconstant(char *s);

%}

letter   [A-Za-z_]
digit    [0-9]
integer  [+-]?[0-9]+
dreal    ([+-]?[0-9]*"."[0-9]*)
ereal    ([+-]?[0-9]*("."[0-9]*)?[eE][+-]?[0-9]+)
real     {dreal}|{ereal}
string   {letter}({letter}|{digit})*

%%
include      BEGIN(incl);

<incl>[ \t]*      /* eat the whitespace */
<incl>[^ \t\n;]+  { /* got the include file name */
        if ( include_stack_ptr >= MAX_INCLUDE_DEPTH )
            {
            fprintf( stderr, "Includes nested too deeply" );
            exit( 1 );
            }

        include_stack[include_stack_ptr++] =
            YY_CURRENT_BUFFER;

        yyin = fopen( yytext, "r" );

        if ( ! yyin ) {
	  printf("ERROR\n");
	  exit (1);
        }

        yy_switch_to_buffer(
            yy_create_buffer( yyin, YY_BUF_SIZE ) );

        BEGIN(INITIAL);
        }

<<EOF>> {
        if ( --include_stack_ptr < 0 )
            {
            yyterminate();
            }

        else
            {
            yy_delete_buffer( YY_CURRENT_BUFFER );
            yy_switch_to_buffer(
                 include_stack[include_stack_ptr] );
            }
        }

[#@](.*)[\n]? {
             extern int lineno;
             lineno++;
             /* we ignore comments */
           }

[ \t]+       /* ignore whitespace */;

[=;,]        return *yytext;

\n         {
             extern int lineno;
             lineno++;
             /* we ignore new line characters */
           }

\{         return LBRACE;

\}         return RBRACE;

input     {
            return INPUT;
          }

{integer} {
            sscanf(yytext, "%d", &yylval.integer);
            return INTEGER;
          }

{real}    {
             sscanf(yytext, "%lf", &yylval.real);
             return REAL;
          }

{string}  {
             int r, pos;

             /* we determine if it is a reserverd word */
	     r = c_symbol_search(yytext, reserved_words,
                         sizeof(reserved_words)/sizeof(reserved_words[0]),
                         &pos);
	     if (r >= 0) {
	       yylval.integer = pos;
	       return IDENTIFIER;
	     }

	     /* otherwise we determine if it is a command */
	     r = c_symbol_search(yytext, commands,
			 sizeof(commands)/sizeof(commands[0]),
			 &pos);
	     if (r >= PRINT) {
	       yylval.integer = r;
	       return r;
	     }

	     /* or a distribution constant */
	     r = c_symbol_search(yytext, distributions,
			 sizeof(distributions)/sizeof(distributions[0]),
			 &pos);
	     if (r >= CONSTANT) {
	       yylval.integer = r;
	       return DISTRIBUTION;
	     }

	     /* or a method */
	     r = c_symbol_search(yytext, methods,
			 sizeof(methods)/sizeof(methods[0]), &pos);
	     if (r >= MCBTR) {
	       yylval.integer = r;
	       return OPTMETHOD;
	     }

             /* or a hessian update method */
             r = c_symbol_search(yytext, hessian_updates,
                         sizeof(hessian_updates)/sizeof(hessian_updates[0]),
			 &pos);
             if (r >= 0) {
               yylval.integer = r;
               return HESSIANUPDATE;
             }

             /* or a random generator */
             r = c_symbol_search(yytext, rangen,
			 sizeof(rangen)/sizeof(rangen[0]),
			 &pos);
	     if (r >= RAN_LECUYER88B) {
	       yylval.integer = r;
	       return RANGENERATOR;
	     }

             /* or a sampling method */
             r = c_symbol_search(yytext, moca_methods,
			 sizeof(moca_methods)/sizeof(moca_methods[0]),
			 &pos);
	     if (r >= 0) {
	       yylval.integer = r;
	       return SAMPLINGMETHOD;
	     }

             /* or a randomization method for Halton sequences */
             r = c_symbol_search(yytext, randomization_halton,
                         sizeof(randomization_halton)/
			 sizeof(randomization_halton[0]), &pos);
             if (r >= MOCA_HALTON_RANDOMIZATION_NONE) {
	       yylval.integer = r;
	       return RANDOMIZATIONMETHOD;
             };

	     /* or a data format */
	     r = c_symbol_search(yytext, formats,
			 sizeof(formats)/sizeof(formats[0]), &pos);
	     if (r > 0) {
	       yylval.integer = r;
	       return FORMAT_DATA;
	     }

             /* or a nonparametric method */
             r = c_symbol_search(yytext, nonparametric_methods,
			 sizeof(nonparametric_methods)/
                         sizeof(nonparametric_methods[0]), &pos);
	     if (r > 0) {
	       yylval.integer = r;
	       return NONPARAMETRIC;
	     }


             /* or a numerical constant */
             r = isconstant(yytext);
             if (r >= 0) {
               yylval.integer = r;
               return INTEGER;
             }

             yylval.string = strdup(yytext);
             return VARIABLE;
          }

\"[^"\n]*["\n"] { /* String constants */
            yylval.string = strdup(yytext+1);
            if (yylval.string[yyleng-2] != '"')
              c_warning("improperly terminated string");
            else
	      yylval.string[yyleng-2] = 0;
	    return STRING;
          }

\[{integer}\] {
                sscanf(yytext, "[%d]", &yylval.integer);
		return POS;
              }

%%

int isconstant(char *s)
{
  if (!strcasecmp(s, "true") || !strcasecmp(s, "new")) return 1;
  if (!strcasecmp(s, "false") || !strcasecmp(s, "old")) return 0;
  return -1;
}

void yyinit(FILE *f) {
  yyin = f;
}

int yywrap(void) {
  return 1;
}
