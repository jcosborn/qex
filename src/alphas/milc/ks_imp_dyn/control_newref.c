/*

Rewrite of "update_onemass.c" that is slightly more clear than 
the legacy code. Also includes 2nd-order Omelyan integrator,
which is better than the leapfrog algorithm that the
legacy code uses. R-algorithm option removed because nobody
uses it anymore...

Author: Curtis Taylor Peterson

*/

#include "ks_imp_includes.h"

void predict_next_xxx(Real *oldtime,Real *newtime,Real *nexttime){
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

int update(){
    // Variables for HMC
    void predict_next_xxx(Real *oldtime,Real *newtime,Real *nexttime);
    double startaction,endaction,d_action(),dH,xrandom;
    int step,iters = 0;
    Real final_rsq;
    Real lmbda = 0.1931833275037836;
    Real cg_time = 0.0;
    Real old_cg_time,next_cg_time = -1.0e6;
    imp_ferm_links_t** fn;

    // Momentum update
    void update_t(int step, Real offset){
        next_cg_time = ((Real)step-offset)*epsilon;
	    predict_next_xxx(&old_cg_time,&cg_time,&next_cg_time);
        restore_fermion_links_from_site(fn_links, MILC_PRECISION);
	    fn = get_fm_links(fn_links);
     	iters += ks_congrad( 
            F_OFFSET(phi), F_OFFSET(xxx), mass, 
			niter, nrestart, rsqmin, MILC_PRECISION, 
			EVEN, &final_rsq, fn[0] 
        );
	    dslash_site( F_OFFSET(xxx), F_OFFSET(xxx), ODD, fn[0]);
	    cg_time = ((Real)step - offset)*epsilon;
        update_h(0.5*epsilon);
    }

    // First smearing
    restore_fermion_links_from_site(fn_links, MILC_PRECISION);
	fn = get_fm_links(fn_links);

    // Heatbath (momentum & fermion)
    ranmom();
    grsource_imp(F_OFFSET(phi), mass, EVEN, fn[0]);

    // Calculate initial action & store field
    iters += ks_congrad( 
        F_OFFSET(phi), F_OFFSET(xxx), mass, 
		niter, nrestart, rsqmin, MILC_PRECISION, 
		EVEN, &final_rsq, fn[0] 
    );
    startaction = d_action();
    gauge_field_copy(F_OFFSET(link[0]), F_OFFSET(old_link[0]));

    // Run HMC iterations: 2nd-order Omelyan
    for(step=1; step <= steps; step++){
        if (step == 1){update_u(lmbda*epsilon);}
        else {update_u(2.0*lmbda*epsilon);}
        update_t(step,0.25);
        update_u((1.0-2.0*lmbda)*epsilon);
        update_t(step,0.5);
        if (step == steps){update_u(lmbda*epsilon);}
    }

    // Calculate final action
    next_cg_time = steps*epsilon;
    predict_next_xxx(&old_cg_time,&cg_time,&next_cg_time);
    restore_fermion_links_from_site(fn_links, MILC_PRECISION);
    fn = get_fm_links(fn_links);
    iters += ks_congrad( 
        F_OFFSET(phi), F_OFFSET(xxx), mass,
		niter, nrestart, rsqmin, MILC_PRECISION, 
		EVEN, &final_rsq, fn[0] 
    );
    cg_time = steps*epsilon;
    endaction = d_action();

    // Metropolis accept/reject
    if(this_node==0){xrandom = myrand(&node_prn);}
    broadcast_float(&xrandom);
    dH = (double)(endaction-startaction);
    if(exp(-dH)<xrandom){
	    if(steps > 0){
	        gauge_field_copy( F_OFFSET(old_link[0]), F_OFFSET(link[0]) );   
            #ifdef FN
	            invalidate_fermion_links(fn_links);
            #endif
	        node0_printf("REJECT: delta S = %e\n",dH);
        }
    }
    else {
        node0_printf("ACCEPT: delta S = %e\n",dH);
        reunitarize_ks();
    }

    // finish up
    if(steps > 0){return (iters/steps);}
    else {return(-99);}
}