/**
 * dsarray.h
 *
 * D Struct
 * Copyright (C) 2000-2002 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __ORATIO_DSARRAY_H__
#define __ORATIO_DSARRAY_H__

#include "dstruct.h"

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Compute the cumulated frequency of the array elements that are
   * inferior to some criterion.
   *
   * @param base is a pointer towards a sorted array of nmemb elements,
   *             each of size "size".
   *        search est une fonction prenant comme premier argument une variable
   *        de même type que les éléments du tableau, et comme second argument
   *        le critère de recherche. search renvoie un entier > 0 si l'élément
   *        à rechercher est inférieur au premier argument (suivant l'ordre
   *        considéré), un entier < 0 s'il est supérieur, et 0 si l'élement
   *        à rechercher est trouvé.
   */
  double ds_array_cum(void *base, size_t nmemb, size_t size, void *criterion,
		      DSCompareFunc search);

#ifdef __cplusplus
}
#endif

#endif  /* __ORATIO_DSARRAY_H__ */
