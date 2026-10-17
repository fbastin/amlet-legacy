/**
 * btr.hpp
 *
 * OPHELIA Library
 * Copyright (C) 2015 Fabian Bastin <fabian.bastin@cirrelt.ca>
 *
 * Licensed under the Open Software License version 2.1
 */

#ifndef __OPHELIA_BTR_HPP__
#define __OPHELIA_BTR_HPP__

class OP_BTR {

  double radius;   // trust-region radius
  
public:

  virtual double ratio();
  
  void update_radius();
};

#endif   /* __OPHELIA_BTR_HPP__ */
