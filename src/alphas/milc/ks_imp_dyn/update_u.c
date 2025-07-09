/*************** update_u.c ************************************/
/* MIMD version 7 */

/* update the link matrices					*
 *  								*
 *  Go to sixth order in the exponential of the momentum		*
 *  matrices, since unitarity seems important.			*
 *  Evaluation is done as:					*
 *	exp(H) * U = ( 1 + H + H^2/2 + H^3/3 ...)*U		*
 *	= U + H*( U + (H/2)*( U + (H/3)*( ... )))		*
 *								*
 */

#include "ks_imp_includes.h"	/* definitions files and prototypes */

#if defined HAVE_QUDA && defined USE_GF_GPU
#include "../include/generic_quda.h"

/* QUDA version */

void update_u(Real eps){

  int i,dir;
  site *s;
  int j;
  int dim[4] = {nx, ny, nz, nt};

#ifdef FN
  invalidate_fermion_links(fn_links);
#endif

  initialize_quda();
  Real *momentum = (Real*)malloc(sites_on_node*4*sizeof(anti_hermitmat));
  Real *gauge = (Real*)malloc(sites_on_node*4*sizeof(su3_matrix));
  // Populate gauge and momentum fields
  FORALLSITES(i,s){
    for(dir=XUP; dir<=TUP; ++dir){
      for(j=0; j<18; ++j){
        gauge[(4*i + dir)*18 + j] = ((Real*)(&(s->link[dir])))[j];
      }
      for(j=0; j<10; ++j){
        momentum[(4*i + dir)*10 + j] = ((Real*)(&(s->mom[dir])))[j];
      }
    } // dir
  }

  qudaUpdateU(MILC_PRECISION, eps, momentum, gauge);

  // Copy updated gauge field back to site structure
  FORALLSITES(i,s){
    for(dir=XUP; dir<=TUP; ++dir){
      for(j=0; j<18; ++j){
        ((Real*)(&(s->link[dir])))[j] = gauge[(4*i + dir)*18 + j];
      }
    }
  }
  free(momentum);
  free(gauge);
  return;
}

#else

/* CPU version */

void update_u( Real eps ){

  register int i,dir;
  register site *s;
  su3_matrix *link,temp1,temp2,htemp;
  register Real t2,t3,t4,t5,t6;
  /**TEMP**
    Real gf_x,gf_av,gf_max;
    int gf_i,gf_j;
   **END TEMP **/

  /**double dtime,dtime2,dclock();**/
  /**dtime = -dclock();**/

  /* Temporary by-hand optimization until pgcc compiler bug is fixed */
  t2 = eps/2.0;
  t3 = eps/3.0;
  t4 = eps/4.0;
  t5 = eps/5.0;
  t6 = eps/6.0;

  /** TEMP **
    gf_av=gf_max=0.0;
   **END TEMP**/
#ifdef FN
  invalidate_fermion_links(fn_links);
#endif

  FORALLSITES(i,s){
    for(dir=XUP; dir <=TUP; dir++){
      uncompress_anti_hermitian( &(s->mom[dir]) , &htemp );
      link = &(s->link[dir]);
      mult_su3_nn(&htemp,link,&temp1);
      /**scalar_mult_add_su3_matrix(link,&temp1,eps/6.0,&temp2);**/
      scalar_mult_add_su3_matrix(link,&temp1,t6,&temp2);
      mult_su3_nn(&htemp,&temp2,&temp1);
      /**scalar_mult_add_su3_matrix(link,&temp1,eps/5.0,&temp2);**/
      scalar_mult_add_su3_matrix(link,&temp1,t5,&temp2);
      mult_su3_nn(&htemp,&temp2,&temp1);
      /**scalar_mult_add_su3_matrix(link,&temp1,eps/4.0,&temp2);**/
      scalar_mult_add_su3_matrix(link,&temp1,t4,&temp2);
      mult_su3_nn(&htemp,&temp2,&temp1);
      /**scalar_mult_add_su3_matrix(link,&temp1,eps/3.0,&temp2);**/
      scalar_mult_add_su3_matrix(link,&temp1,t3,&temp2);
      mult_su3_nn(&htemp,&temp2,&temp1);
      /**scalar_mult_add_su3_matrix(link,&temp1,eps/2.0,&temp2);**/
      scalar_mult_add_su3_matrix(link,&temp1,t2,&temp2);
      mult_su3_nn(&htemp,&temp2,&temp1);
      scalar_mult_add_su3_matrix(link,&temp1,eps    ,&temp2); 
      su3mat_copy(&temp2,link);
    }
  }

  /**dtime += dclock();
    node0_printf("LINK_UPDATE: time = %e  mflops = %e\n",
    dtime, (double)(5616.0*volume/(1.0e6*dtime*numnodes())) );**/
} /* update_u */

#endif
