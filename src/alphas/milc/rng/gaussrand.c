/*****************  gaussrand.c  (in su3.a) *****************************
 **
 *  Real gaussian_ran_no( double_prn *prn_pt )*
 *  Gaussian distributed random number*
 *  Probability distribution exp( -x*x ), so < x^2 > = 1/2*
 *  This requires a random number generator named "myrand()", returning*
 *  a Real uniformly distributed between zero and one. The argument of*
 *  this routine is a pointer to be passed to myrand(). *
 *  Y. Shamir 9/8/13 - use prn_pt to toggle -- allows
 *  layout-independent results under SITE_PRN when generating an odd
 *  number of values on a site before going to a new site.
*   DT 12/15/15 save state in double_prn structure.  This should be thread-safe
*   when used in a loop over lattice sites (though not if you loop over lots
*   of random numbers in the same site using OpenMP)
 */

#include "../include/config.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/su3.h"
#include "../include/random.h"

//Real gaussian_rand_no( double_prn *prn_pt ){
//  Real fac,r,v1,v2;
//
//  if  (prn_pt->gauss_toggle==0) {
//    do {
//      v1=2.0*myrand(prn_pt)-1.0;
//      v2=2.0*myrand(prn_pt)-1.0;
//      r=v1*v1+v2*v2;
//    } while (r >= 1.0);
//    fac=sqrt( -log((double)r)/(double)r);
//    prn_pt->gauss_gset=v1*fac;
//    prn_pt->gauss_toggle=1;;
//    return v2*fac;
//  } else {
//    prn_pt->gauss_toggle=0;
//    return prn_pt->gauss_gset;
//  }
//}

Real gaussian_rand_no( double_prn *prn_pt ){
  static double_prn *ptoggle=NULL;
  static Real gset;
  Real fac,r,v1,v2;

  if  (ptoggle != prn_pt) {
    do {
      v1=2.0*myrand(prn_pt)-1.0;
      v2=2.0*myrand(prn_pt)-1.0;
      r=v1*v1+v2*v2;
    } while (r >= 1.0);
    fac=sqrt( -log((double)r)/(double)r);
    gset=v1*fac;
    ptoggle=prn_pt;
    return v2*fac;
  } else {
    ptoggle=NULL;
    return gset;
  }
}

complex complex_gaussian_rand_no( double_prn *prn_pt ){
  complex result;

  result.real = gaussian_rand_no( prn_pt );
  result.imag = gaussian_rand_no( prn_pt );
  return result;
}

