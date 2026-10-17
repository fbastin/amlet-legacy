/**
 * btr.hpp
 *
 * OPHELIA Library
 * Copyright (C) 2015 Fabian Bastin <fabian.bastin@gmail.com>
 *
 * Licensed under the MIT License
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
