/*************************** ranmom.c *******************************/
/* MIMD version 7 */

/* Produce Gaussian random momenta for the gauge fields. */

#include "generic_includes.h"
#include <defines.h>                 /* For SITERAND */
#include "../include/generic_quda.h"

void ranmom(){
register int i,dir;
register site *s;
    FORALLSITES(i,s){
	for(dir=XUP;dir<=TUP;dir++){
#ifdef SCHROED_FUN
	    if(dir==TUP || s->t>0){
#endif
#ifdef SITERAND
		random_anti_hermitian( (anti_hermitmat *)&(s->mom[dir]),
		    &(s->site_prn) );
#else
		random_anti_hermitian( (anti_hermitmat *)&(s->mom[dir]),
		    &node_prn );
#endif
#ifdef SCHROED_FUN
	    }
	    else{
		s->mom[dir].m00im = 0.0;
		s->mom[dir].m11im = 0.0;
		s->mom[dir].m22im = 0.0;
		s->mom[dir].m01.real = 0.0;
		s->mom[dir].m01.imag = 0.0;
		s->mom[dir].m02.real = 0.0;
		s->mom[dir].m02.imag = 0.0;
		s->mom[dir].m12.real = 0.0;
		s->mom[dir].m12.imag = 0.0;
	    }
#endif
	}
    }

#if defined(HAVE_QUDA) && defined(USE_GF_GPU) && defined(USE_FF_GPU)
    // if gauge force and fermion force are offloaded then we can keep the momentum resident
    QudaMILCSiteArg_t arg = newQudaMILCSiteArg();
    qudaMomLoad(MILC_PRECISION, &arg);
#endif

}

