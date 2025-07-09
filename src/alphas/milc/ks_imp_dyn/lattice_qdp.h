#ifndef _LATTICE_QDP_H
#define _LATTICE_QDP_H

#ifdef HAVE_QDP
#include <qdp.h>

#ifdef CONTROL
#define EXTERN 
#else
#define EXTERN extern
#endif

EXTERN QDP_Shift neighbor3[4];
EXTERN QDP_Shift shiftdirs[8];
EXTERN QDP_ShiftDir shiftfwd[8], shiftbck[8];

#endif
#endif /* _LATTICE_QDP_H */
