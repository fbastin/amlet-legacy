/**
 * exception.h
 *
 * C Tools
 * Copyright (C) 2000-2002 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

/**
 * Allows exceptions in C.
 *
 * The follwoing code has been inspired by the paper of Jakub Zimmerman,
 * published in Login (68), December 1999.
 */

#ifndef __CTOOLS_EXCEPTION_H__
#define __CTOOLS_EXCEPTION_H__

#include <setjmp.h>

typedef struct except_ctx {
  jmp_buf ctx;
  struct except_ctx *super;
} except;

extern except *_exceptions;

#define TRY { \
         int exc; \
         push_exc(); \
         exc=setjmp(_exceptions->ctx); \
         switch(exc) { \
	 case 0: 

#define CATCH(X) \
         break; \
	 case (X): \
	 exc=0;

#define ENDTRY \
         break; \
         default :\
         } \
         pop_exc(exc); \
         }

#define THROW(X) longjmp(_exceptions->ctx, X)

#endif  /* __CTOOLS_EXCEPTION_H__ */
