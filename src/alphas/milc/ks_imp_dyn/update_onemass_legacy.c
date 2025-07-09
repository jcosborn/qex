/********** update_onemass.c ***************************************/
/* MIMD version 7 */

/*
 Update lattice.
 Improved method for 1-4 flavors:
	update U by (epsilon/2)*(1-Nf/4)
	compute PHI
	update U to epsilon/2
	compute X
	update H, full step
	update U to next time needed

 This routine does not refresh the antihermitian momenta.
 This routine begins at "integral" time, with H and U evaluated
 at same time.
*/
#include "ks_imp_includes.h"	/* definitions files and prototypes */

int update()  {
int step, iters=0;
Real final_rsq;
void predict_next_xxx(Real *oldtime,Real *newtime,Real *nexttime);
Real cg_time;	/* simulation time for last two CG's */
#ifdef HMC_ALGORITHM
Real old_cg_time,next_cg_time;	/* simulation time for last two CG's */
double startaction,endaction,d_action();
Real xrandom;
Real lmbda = 0.1931833275037836;
#endif
  imp_ferm_links_t** fn;

    /* refresh the momenta */
    ranmom();

/*DEBUG*/
/**grsource_imp(F_OFFSET(phi), mass, EVENANDODD);
ks_congrad(F_OFFSET(phi),F_OFFSET(xxx),mass,niter,nrestart,rsqmin,MILC_PRECISION,EVENANDODD,&final_rsq);
checkmul();**/
/*ENDDEBUG*/

    /* do "steps" microcanonical steps  */
    for(step=1; step <= steps; step++){
 
#ifdef PHI_ALGORITHM
        /* generate a pseudofermion configuration only at start*/
      if(step==1){
	restore_fermion_links_from_site(fn_links, MILC_PRECISION);
	fn = get_fm_links(fn_links);
	grsource_imp(F_OFFSET(phi), mass, EVEN, fn[0]); 
	old_cg_time = cg_time = -1.0e6;
      }

#ifdef HMC_ALGORITHM
        /* find action */
        /* do conjugate gradient to get (Madj M)inverse * phi */
        if(step==1){
            /* do conjugate gradient to get (Madj M)inverse * phi */
	  restore_fermion_links_from_site(fn_links, MILC_PRECISION);
	  fn = get_fm_links(fn_links);
	    iters += ks_congrad( F_OFFSET(phi), F_OFFSET(xxx), mass, 
				 niter, nrestart, rsqmin, MILC_PRECISION, 
				 EVEN, &final_rsq, fn[0] );
	    cg_time = 0.0;

     	    startaction=d_action();
            /* copy link field to old_link */
	    gauge_field_copy( F_OFFSET(link[0]), F_OFFSET(old_link[0]));
        }
#endif

	/* update U's to middle of interval */
     	//update_u(0.5*epsilon);
    update_u(lmbda*epsilon);

	/* save conjugate gradient solution, predict next one */
	next_cg_time = ((Real)step-0.5)*epsilon;
	predict_next_xxx(&old_cg_time,&cg_time,&next_cg_time);

#else /* "R" algorithm */
       	/* first update the U's to special time interval */
       	update_u(epsilon*(0.5-nflavors/8.0));

        /* generate a pseudofermion configuration */
	restore_fermion_links_from_site(fn_links, MILC_PRECISION);
	fn = get_fm_links(fn_links);
     	grsource_imp(F_OFFSET(phi), mass, EVEN, fn[0]); 
	cg_time = -1.0e6;

	/* update U's to middle of interval */
     	update_u(epsilon*nflavors/8.0);
#endif

        /* do conjugate gradient to get (Madj M)inverse * phi */
	restore_fermion_links_from_site(fn_links, MILC_PRECISION);
	fn = get_fm_links(fn_links);
     	iters += ks_congrad( F_OFFSET(phi), F_OFFSET(xxx), mass, 
			     niter, nrestart, rsqmin, MILC_PRECISION, 
			     EVEN, &final_rsq, fn[0] );
	dslash_site( F_OFFSET(xxx), F_OFFSET(xxx), ODD, fn[0]);
	cg_time = ((Real)step - 0.333)*epsilon;
	/* now update H by full time interval */
    	update_h(0.5*epsilon);

    	/* update U's by half time step to get to even time */
    	update_u((1.0-2.0*lmbda)*epsilon); //update_u(epsilon*0.5);

        /* do conjugate gradient to get (Madj M)inverse * phi */
	restore_fermion_links_from_site(fn_links, MILC_PRECISION);
	fn = get_fm_links(fn_links);
     	iters += ks_congrad( F_OFFSET(phi), F_OFFSET(xxx), mass, 
			     niter, nrestart, rsqmin, MILC_PRECISION, 
			     EVEN, &final_rsq, fn[0] );
	dslash_site( F_OFFSET(xxx), F_OFFSET(xxx), ODD, fn[0]);
	cg_time = ((Real)step - 0.333)*epsilon;
	/* now update H by full time interval */
    	update_h(0.5*epsilon);

    update_u(lmbda*epsilon);

        /* reunitarize the gauge field */
        reunitarize_ks();


    }	/* end loop over microcanonical steps */

#ifdef HMC_ALGORITHM
    /* find action */
    /* do conjugate gradient to get (Madj M)inverse * phi */
    next_cg_time = steps*epsilon;
    predict_next_xxx(&old_cg_time,&cg_time,&next_cg_time);
    restore_fermion_links_from_site(fn_links, MILC_PRECISION);
    fn = get_fm_links(fn_links);
    iters += ks_congrad( F_OFFSET(phi), F_OFFSET(xxx), mass,
			 niter, nrestart, rsqmin, MILC_PRECISION, 
			 EVEN, &final_rsq, fn[0] );
    cg_time = steps*epsilon;
    endaction=d_action();
    /* decide whether to accept, if not, copy old link field back */
    /* careful - must generate only one random number for whole lattice */
    if(this_node==0)xrandom = myrand(&node_prn);
    broadcast_float(&xrandom);
    if( exp( (double)(startaction-endaction) ) < xrandom ){
	if(steps > 0)
	    gauge_field_copy( F_OFFSET(old_link[0]), F_OFFSET(link[0]) );
#ifdef FN
	invalidate_fermion_links(fn_links);
#endif
	node0_printf("REJECT: delta S = %e\n", (double)(endaction-startaction));
    }
    else {
	node0_printf("ACCEPT: delta S = %e\n", (double)(endaction-startaction));
    }
#endif

    if(steps > 0)return (iters/steps);
    else return(-99);
}

#ifdef HMC_ALGORITHM
/* use linear extrapolation to predict next conjugate gradient solution */
/* only need even sites */
void predict_next_xxx(Real *oldtime,Real *newtime,Real *nexttime)
{
register int i;
register site *s;
register Real x;
su3_vector tvec;
    if( *newtime != *oldtime ) x = (*nexttime-*newtime)/(*newtime-*oldtime);
    else x = 0.0;
    if( *oldtime < 0.0 ){
        FOREVENSITES(i,s){
	    s->old_xxx = s->xxx;
        }
    }
    else  {
        FOREVENSITES(i,s){
            sub_su3_vector( &(s->xxx), &(s->old_xxx), &tvec);
	    s->old_xxx = s->xxx;
            scalar_mult_add_su3_vector( &(s->xxx), &tvec,x, &(s->xxx) );
        }
    }
    *oldtime = *newtime;
    *newtime = *nexttime;
}
#endif
