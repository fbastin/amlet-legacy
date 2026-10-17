/**
 * nt_coordinates.h
 *
 * ORATIO Library
 * Copyright (C) 2011 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under BSD License.
 */

#ifndef __ORATIO_NTOOLS_COORDINATES_H__
#define __ORATIO_NTOOLS_COORDINATES_H__

#ifdef __cplusplus
extern "C" {
#endif
	
	/**
	 * Convert cartesian coordinates to 2D polar coordinates.
	 *
	 * @param x is a 2 dimensional vector that contains the
	 *   Cartesian coordinates.
	 * @param r is on output the radial distance from the origin.
	 * @param theta is on output the counterclockwise angle from the
	 *   x-axis.
	 */
	void nt_cartesian2polar(double *x, double r, double theta);
	
  /**
   * Convert 2D polar coordinates in cartesian coordinates.
   *
   * @param r is the radial distance from the origin.
   * @param theta is the counterclockwise angle from the x-axis.
   * @param x is a 2 dimensional vector that will contain the
   *   Cartesian coordinates on output.
   */
  void nt_polar2cartesian(double r, double theta, double *x);
	
#ifdef __cplusplus
}
#endif

#endif /* __ORATIO_NT_COORDINATES_H__ */
