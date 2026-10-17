/**
 * machar_test.c
 *
 * ORATIO Library
 * Copyright (C) 1999-2010 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#include <stdio.h>

#include <oratio/ntools.h>

/**
 * This routine prints hardware-determined double-precision
 * machine constants obtained from nt_machar.  Dmachar is a C
 * translation of the Fortran routine MACHAR from W. J. Cody,
 * "MACHAR: A subroutine to dynamically determine machine
 * parameters".  TOMS (14), 1988.
 * 
 * Descriptions of the machine constants are given in the
 * prologue comments in rmachar.
 *  
 * Original driver: Richard Bartels, October 16, 1985
 * 
 * Modified by: W. J. Cody
 * July 26, 1988
 **/
int nt_machar_test()
{
  int ibeta,iexp,irnd,it,machep,maxexp,minexp,negep,ngrd;
  int i ;
  double eps,epsneg,xmax,xmin;
  union wjc{
    long int jj[2];
    double xbig;
  } uval;

  nt_machar(&ibeta,&it,&irnd,&ngrd,&machep,&negep,
	    &iexp,&minexp,&maxexp,&eps,&epsneg,&xmin,&xmax);
  printf("Double precision MACHAR constants\n");
  printf("ibeta  = %d\n",ibeta);
  printf("it     = %d\n",it);
  printf("irnd   = %d\n",irnd);
  printf("ngrd   = %d\n",ngrd);
  printf("machep = %d\n",machep);
  printf("negep  = %d\n",negep);
  printf("iexp   = %d\n",iexp);
  printf("minexp = %d\n",minexp);
  printf("maxexp = %d\n",maxexp);

#define DISPLAY(s,x) { \
	uval.xbig = x ; \
	printf(s); \
	printf(" %24.16e ",(double) x) ; \
	for(i=0;i<2;i++) printf(" %9lX ",uval.jj[i]) ; \
	printf("\n"); \
	}
			
  DISPLAY("eps   ",eps);
  DISPLAY("epsneg",epsneg);
  DISPLAY("xmin  ",xmin);
  DISPLAY("xmax  ",xmax);

  return 0;
}
