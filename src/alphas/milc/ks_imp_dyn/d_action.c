/*************** d_action.c ****************************************/
/* MIMD version 7 */

/* Measure total action, as needed by the hybrid Monte Carlo
   algorithm.  When this routine is called the conjugate gradient
   should already have been run on the even sites, so that in the case
   of two masses, the vectors xxx1 and xxx2 contain (M_adjoint*M)^(-1)
   * phi1 and phi2 and in the case of one mass the vector xxx contains
   (M_adjoint*M)^(-1) * phi.
*/

#include "ks_imp_includes.h"	/* definitions files and prototypes */
Real ahmat_mag_sq(anti_hermitmat *pt);

/*DEBUG*/
double old_g, old_h, old_f, old_a;
/*ENDDEBUG*/

double d_action(){
  double hmom_action(),fermion_action();
  double ssplaq,stplaq,g_action,h_action,f_action;

  h_action = hmom_action();

  d_plaquette(&ssplaq,&stplaq);
  ssplaq *= -1.0; stplaq *= -1.0;
  g_action = -beta*volume*(ssplaq+stplaq);
  node0_printf("PLAQUETTE ACTION: %e\n",g_action);
  rephase(OFF);
  g_action = (beta/3.0)*imp_gauge_action();
  rephase(ON);

  f_action = fermion_action();

  node0_printf("ACTION: g,h,f = %0.17lf  %0.17lf  %0.17lf  %0.17lf\n",
  g_action, h_action, f_action, g_action+h_action+f_action );

  node0_printf("DG = %0.17lf, DH = %0.17lf, DF = %0.17lf, D = %0.17lf\n",
  g_action-old_g, h_action-old_h, f_action-old_f,
  g_action+h_action+f_action-old_a);
  old_g=g_action; old_h=h_action; old_f=f_action;
  old_a=g_action+h_action+f_action;

  return(g_action+h_action+f_action);
}

/* fermion contribution to the action */
double fermion_action() {
  register int i;
  register site *s;
  register complex cc;
  register complex hcc;
  double sum = 0.0;

  /*
  #ifdef HASENBUSCH
      dslash_site(F_OFFSET(hxxx), F_OFFSET(hxxxr), EVEN, fn[0]);
  #endif
  dslash_site(F_OFFSET(xxx), F_OFFSET(xxxr), EVEN, fn[0]);
  */

  FORALLSITES(i,s){
    #ifdef ONEMASS
      #ifdef HASENBUSCH
        //hcc = su3_dot( &(s->hphi), &(s->hxxx) );
        //sum += 4.0*(hmass*hmass - mass*mass)*(double)hcc.real; // det(M)/det(H)
        hcc = su3_dot(&(s->hpsi), &(s->hpsi));
        sum += (double)hcc.real;
      #endif
      //cc = su3_dot( &(s->phi), &(s->xxx) );
      //sum += (double)cc.real; // det(H)
      cc = su3_dot(&(s->psi), &(s->psi));
      sum += (double)cc.real; // det(H)
    #else
      cc = su3_dot( &(s->phi1), &(s->xxx1) );
      sum += (double)cc.real;
      cc = su3_dot( &(s->phi2), &(s->xxx2) );
      sum += (double)cc.real;
    #endif
  }
  g_doublesum(&sum);
  return (sum);
}

/* gauge momentum contribution to the action */
double hmom_action() {
  register int i,dir;
  register site *s;
  double sum;

  sum = 0.0;
  FORALLSITES(i,s){
	for(dir=XUP;dir<=TUP;dir++){
        sum += (double)ahmat_mag_sq( &(s->mom[dir]) ) - 4.0;
	}
  }
  g_doublesum( &sum );
  return(sum);
}

/* magnitude squared of an antihermition matrix */
Real ahmat_mag_sq(anti_hermitmat *pt){
  register Real x,sum;
  x = pt->m00im; sum  = 0.5*x*x;
  x = pt->m11im; sum += 0.5*x*x;
  x = pt->m22im; sum += 0.5*x*x;
  x = pt->m01.real; sum += x*x;
  x = pt->m01.imag; sum += x*x;
  x = pt->m02.real; sum += x*x;
  x = pt->m02.imag; sum += x*x;
  x = pt->m12.real; sum += x*x;
  x = pt->m12.imag; sum += x*x;
  return(sum);
}

