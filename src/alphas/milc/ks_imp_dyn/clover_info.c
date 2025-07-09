/*********************** clover_info.c *************************/
/* MIMD version 7 */

/* For clover_invert */

/* Application-dependent routine for writing gauge info file
   called from one of the output routines in io_prop_w.c */

/* This file is an ASCII companion to the gauge configuration file
   and contains information about the action used to generate it.
   This information is consistently written in the pattern

       keyword  value

   To maintain a semblance of consistency, the possible keywords are
   listed in io_wprop.h.  Add more as the need arises, but be sure
   to notify the rest of the collaboration.

   */

/* build_w_prop_hdr        Fills in the spin table of contents in the header
                           structure
   write_appl_w_prop_info  Writes supplementary information to the info file */

#include "ks_imp_includes.h"
#include <string.h>

/*---------------------------------------------------------------------------*/

void write_appl_w_prop_info(FILE *fp)
{

  /* Note: this application does not write propagator files! */
}

#define INFOSTRING_MAX 2048

char *create_w_QCDML(){
  char *info = NULL;
  return info;
}

void free_w_QCDML(char *info){
  if(info != NULL)free(info);
}

