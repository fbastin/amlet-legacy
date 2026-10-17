/**
 * info.h
 *
 * AMLET (Another Mixed Logit Estimation Tool)
 * Copyright (C) 2002-2009 Fabian Bastin <bastin@iro.umontreal.ca>
 *
 * Licensed under the MIT License
 */

#ifndef __AMLET_INFO_H__
#define __AMLET_INFO_H__

#include <oratio/ctools.h>

#include "amlet_driver.h"

#define COPYRIGHT   "Copyright (C) 2002-2015"
#define AUTHORS     "Fabian Bastin and contributors"
#define EMAIL       "fabian.bastin@gmail.com"
#define COMMAND     "amlet"
#define PROGRAM     "AMLET"

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Error identifiers.
   */

  enum {
    AMLET_NORMAL_EXIT = 0,
    AMLET_HELP_EXIT,
    AMLET_VERSION_EXIT,
    AMLET_ERROR_OPTIONS,
    AMLET_ERROR_CHOICES,
    AMLET_ERROR_PARSING,
    AMLET_ERROR_DRAWS_FILE_OPENING,
    AMLET_ERROR_DRAWS_FILE_READING,
    AMLET_ERROR_REPETITIONS,
    AMLET_ERROR_REPETITIONS_FILE_OPENING,
    AMLET_ERROR_REPETITIONS_FILE_READING,
    AMLET_ERROR_PEOPLE,
    AMLET_ERROR_STARTING_POINT,
    /* Data file */
    AMLET_ERROR_UNKNOWN_FORMAT,
    AMLET_ERROR_FILE_OPEN,
    AMLET_ERROR_DATA_READING,
    AMLET_ERROR_OPTIMISATION,
    /* XML */
    AMLET_ERROR_XML_PARSING
  };

  static const c_symbol amlet_error_msgs[] = {
    {"normal termination", AMLET_NORMAL_EXIT},
    {"", AMLET_HELP_EXIT},
    {"", AMLET_VERSION_EXIT},
    {"error while parsing amlet options", AMLET_ERROR_OPTIONS},
    {"", AMLET_ERROR_CHOICES},
    {"", AMLET_ERROR_PARSING},
    {"unable to open draws file", AMLET_ERROR_DRAWS_FILE_OPENING},
    {"error during draws reading", AMLET_ERROR_DRAWS_FILE_READING},
    {"error with the repetitions number",
     AMLET_ERROR_REPETITIONS},
    {"error while opening repetitions file",
     AMLET_ERROR_REPETITIONS_FILE_OPENING},
    {"error while reading repetitions file",
     AMLET_ERROR_REPETITIONS_FILE_READING},
    {"", AMLET_ERROR_PEOPLE},
    {"", AMLET_ERROR_STARTING_POINT},
    /* Data file */
    {"", AMLET_ERROR_UNKNOWN_FORMAT},
    {"", AMLET_ERROR_FILE_OPEN},
    {"error in data reading", AMLET_ERROR_DATA_READING},
    {"", AMLET_ERROR_OPTIMISATION},
    /* XML */
    {"", AMLET_ERROR_XML_PARSING}
  };

  /**
   * @param info is the numerical identifier of the message to print.
   */
  char *amlet_info_message(int info);

  /**
   * Information Routines.
   */

  /**
   * Print a message describing the basic usage of the program.
   *
   * @param s is the string to display.
   */
  void amlet_help(char *s);

#ifdef __cplusplus
}
#endif

#endif  /* __AMLET_INFO_H__ */
