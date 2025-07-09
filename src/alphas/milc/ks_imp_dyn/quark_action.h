#ifndef _HISQ_U3_ACTION_H
#define _HISQ_U3_ACTION_H

#include "../include/umethod.h"
#include "../include/dirs.h"
#include "../generic_ks/imp_actions/imp_action_types.h"
#define FERM_ACTION HISQ

    /* Specify paths in orientation in which they appear in the
       forward part of the x component of dslash().  Rotations and
       reflections will be automatically included. Be careful
       about signs of coefficients.  See long comment at bottom
       of quark_stuff.c. */

//#define HISQ_NAIK_ADJUSTABLE // allow for adjustable epsilon in Naik term
//#define HISQ_NAIK_2ND_ORDER (-0.675)

#define QUARK_ACTION_DESCRIPTION "\"HISQ action version 1\""

// Smearing for first level
// This is Fat7
//#define ASQ_OPTIMIZED_FORCE_1
#define ASQ_OPTIMIZED_FATTENING
#define QUARK_ACTION_DESCRIPTION_1 "\"Fat 7 (level 1)\""

#define MAX_NUM 688  // should be obsolete, for now max of MAX_NUM_[12]
#define MAX_LENGTH 7	// Maximum length of path in any path table
#define MAX_BASIC_PATHS 6  // Max. no. of basic paths in any path table
#define NUM_BASIC_PATHS_1 4
#define MAX_NUM_1 632
#ifdef IMP_QUARK_ACTION_DEFINE_PATH_TABLES
    static int path_ind_1[NUM_BASIC_PATHS_1][MAX_LENGTH] = {
    { XUP, NODIR, NODIR, NODIR, NODIR, NODIR, NODIR },  /* One Link */
    { YUP, XUP, YDOWN, NODIR, NODIR, NODIR, NODIR },    /* Staple */
    { YUP, ZUP, XUP, ZDOWN, YDOWN, NODIR, NODIR },      /* 5-link for flavor sym. */
    { YUP, ZUP, TUP, XUP, TDOWN, ZDOWN, YDOWN}, /* 7-link for flavor sym. */
    };
    static int quark_action_npaths_1 = NUM_BASIC_PATHS_1 ;
    static int path_length_in_1[NUM_BASIC_PATHS_1] = {1,3,5,7};
    static Real path_coeff_1[NUM_BASIC_PATHS_1] = {
       ( 1.0/8.0),        /* one link */
       (-1.0/8.0)*0.5,              /* simple staple */
       ( 1.0/8.0)*0.25*0.5,         /* displace link in two directions */
       (-1.0/8.0)*0.125*(1.0/6.0),  /* displace link in three directions */
    };
#endif

//#define UNITARIZATION_METHOD UNITARIZE_NONE
//#define UNITARIZATION_METHOD UNITARIZE_ROOT
//#define UNITARIZATION_METHOD UNITARIZE_RATIONAL
#define UNITARIZATION_METHOD UNITARIZE_ANALYTIC

//#define UNITARIZATION_GROUP UNITARIZE_SU3
#define UNITARIZATION_GROUP UNITARIZE_U3

// Smearing for second level
//#define ASQ_OPTIMIZED_FATTENING_2
//#define ASQ_OPTIMIZED_FORCE_2
#define QUARK_ACTION_DESCRIPTION_2 "\"Fat7 + 2xLepage\""

#ifdef IMP_QUARK_ACTION_DEFINE_PATH_TABLES
#define NUM_BASIC_PATHS_2 6
#define MAX_NUM_2 688
    static int path_ind_2[NUM_BASIC_PATHS_2][MAX_LENGTH] = {
    { XUP, NODIR, NODIR, NODIR, NODIR, NODIR, NODIR },  /* One Link */
    { XUP, XUP, XUP, NODIR, NODIR, NODIR, NODIR },      /* Naik */
    { YUP, XUP, YDOWN, NODIR, NODIR, NODIR, NODIR },    /* Staple */
    { YUP, ZUP, XUP, ZDOWN, YDOWN, NODIR, NODIR },      /* 5-link for flavor sym. */
    { YUP, ZUP, TUP, XUP, TDOWN, ZDOWN, YDOWN}, /* 7-link for flavor sym. */
    { YUP, YUP, XUP, YDOWN, YDOWN, NODIR, NODIR },      /* 5-link compensation    */
    };
    static int path_length_in_2[NUM_BASIC_PATHS_2] = {1,3,3,5,7,5};
    static int quark_action_npaths_2 = NUM_BASIC_PATHS_2 ;
    static Real path_coeff_2[NUM_BASIC_PATHS_2] = {
       (( 1.0/8.0)+(2.0*6.0/16.0)+(1.0/8.0)),        /* one link */
            /*One link is 1/8 as in fat7 + 2*3/8 for Lepage + 1/8 for Naik */
       (-1.0/24.0),                 /* Naik */
       (-1.0/8.0)*0.5,              /* simple staple */
       ( 1.0/8.0)*0.25*0.5,         /* displace link in two directions */
       (-1.0/8.0)*0.125*(1.0/6.0),  /* displace link in three directions */
       (-2.0/16 ),                  /* Correct O(a^2) errors, 2X as much as Asqtad  */
    };
#define INDEX_ONELINK 0
#define INDEX_NAIK 1
#endif


// Smearing for second level -- only difference in 1-link and Naik
#define NUM_BASIC_PATHS_3 2
#define MAX_NUM_3 16
#define QUARK_ACTION_DESCRIPTION_3 "\"1-link + Naik\""
#ifdef IMP_QUARK_ACTION_DEFINE_PATH_TABLES
    static int path_ind_3[NUM_BASIC_PATHS_3][MAX_LENGTH] = {
    { XUP, NODIR, NODIR, NODIR, NODIR, NODIR, NODIR },  /* One Link */
    { XUP, XUP, XUP, NODIR, NODIR, NODIR, NODIR },      /* Naik */
    };
    static int path_length_in_3[NUM_BASIC_PATHS_3] = {1,3};
    static int quark_action_npaths_3 = NUM_BASIC_PATHS_3 ;
    static Real path_coeff_3[NUM_BASIC_PATHS_3] = {
       1.0/8.0,        /* one link */
       -1.0/24.0,                 /* Naik */
    };
#endif /* IMP_QUARK_ACTION_INFO_ONLY */


#endif // _HISQ_U3_ACTION_H
