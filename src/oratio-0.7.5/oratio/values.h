/**
 * values.h
 *
 * C Tools
 * Copyright (C) 2005 Fabian Bastin <slashbin@namurlug.org>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __CTOOLS_VALUES_H__
#define __CTOOLS_VALUES_H__

/**
 * Define some system constant values that using the names defined in
 * the GNU C library (file values.h).
 */

#include <limits.h>

#ifndef _VALUES_H

#define CTOOLS_BITS(type) (sizeof (type)*CHAR_BIT)

#ifndef CHARBITS
#define CHARBITS	CTOOLS_BITS (char)
#endif

#ifndef SHORTBITS
#define SHORTBITS	CTOOLS_BITS (short int)
#endif

#ifndef INTBITS
#define INTBITS		CTOOLS_BITS (int)
#endif

#ifndef LONGBITS
#define LONGBITS	CTOOLS_BITS (long int)
#endif

#ifndef PTRBITS
#define PTRBITS		CTOOLS_BITS (char *)
#endif

#ifndef DOUBLEBITS
#define DOUBLEBITS	CTOOLS_BITS (double)
#endif

#ifndef FLOATBITS
#define FLOATBITS	CTOOLS_BITS (float)
#endif

#ifndef DBL_MAX
#include <float.h>
#endif

#endif /* _VALUES_H */

#endif  /* __CTOOLS_VALUES_H__ */
