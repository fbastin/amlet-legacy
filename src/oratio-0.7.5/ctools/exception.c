/**
 *
 *   C Tools
 *   Copyright (C) 2000-2002 Fabian Bastin <slashbin@namurlug.org>
 *
 *   Licensed under the MIT License
 */

#include <stdio.h>
#include <stdlib.h>
#include <oratio/ctools.h>

except *_exceptions = NULL;

void push_exc() {
  except *layer;
  layer = _malloc(except);
  layer->super = _exceptions;
  _exceptions = layer;
}

void pop_exc(int ex) {
  except *super;
  super = _exceptions->super;
  free(_exceptions);
  _exceptions = super;
  if (ex) {
    if (_exceptions)
      longjmp(_exceptions->ctx, ex);
    else {
      fprintf(stderr, "Unrecoverable exception: %d\n", ex);
      exit(-1);
    }
  }
}
