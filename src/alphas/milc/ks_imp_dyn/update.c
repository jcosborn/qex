/********** update.c ****************************************************/
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

#if 1
#ifdef HAVE_QIO
#include <qio.h>
#endif
#endif

int update()  {
  int step, iters=0;
  int n;
  Real final_rsq;
#ifdef HMC_ALGORITHM
  double startaction,endaction,d_action();
  Real xrandom;
#endif
  imp_ferm_links_t** fn;

    /* refresh the momenta */
    ranmom();

    /* In this application, the number of naik terms is 1 or 2 only */
    n = fermion_links_get_n_naiks(fn_links);

    /* do "steps" microcanonical steps"  */
    for(step=1; step <= steps; step++){
 
#ifdef PHI_ALGORITHM
        /* generate a pseudofermion configuration only at start*/
	/* also clear xxx, since zero is our best guess for the solution
	   with a new random phi field. */
     	if(step==1){
	  restore_fermion_links_from_site(fn_links, MILC_PRECISION);
	  fn = get_fm_links(fn_links);
	  clear_latvec( F_OFFSET(xxx1), EVENANDODD );
	  grsource_imp( F_OFFSET(phi1), mass1, EVEN, fn[0]);
	  clear_latvec( F_OFFSET(xxx2), EVENANDODD );
	  grsource_imp( F_OFFSET(phi2), mass2, EVEN, fn[n-1]);
	}

#ifdef HMC_ALGORITHM
        /* find action */
        /* do conjugate gradient to get (Madj M)inverse * phi */
        if(step==1){
            /* do conjugate gradient to get (Madj M)inverse * phi */
	  restore_fermion_links_from_site(fn_links, MILC_PRECISION);
	  fn = get_fm_links(fn_links);
	    iters += ks_congrad( F_OFFSET(phi1), F_OFFSET(xxx1), mass1,
				 niter, nrestart, rsqmin, MILC_PRECISION, EVEN, 
				 &final_rsq, fn[0]);
	  restore_fermion_links_from_site(fn_links, MILC_PRECISION);
	  fn = get_fm_links(fn_links);
	    iters += ks_congrad( F_OFFSET(phi2), F_OFFSET(xxx2), mass2,
				 niter, nrestart, rsqmin, MILC_PRECISION, EVEN, 
				 &final_rsq, fn[n-1]);

     	    startaction=d_action();
            /* copy link field to old_link */
	    gauge_field_copy( F_OFFSET(link[0]), F_OFFSET(old_link[0]));
        }
#endif

	/* update U's to middle of interval */
     	update_u(0.5*epsilon);

#else /* "R" algorithm */
       	/* first update the U's to special time interval */
        /* and generate a pseudofermion configuration */
	/* probably makes most sense if nflavors1 >= nflavors2 */

       	update_u(epsilon*(0.5-nflavors1/8.0));
	clear_latvec( F_OFFSET(xxx1), EVENANDODD );
	restore_fermion_links_from_site(fn_links, MILC_PRECISION);
	fn = get_fm_links(fn_links);
     	grsource_imp( F_OFFSET(phi1), mass1, EVEN, fn[0]);

       	update_u(epsilon*((nflavors1-nflavors2)/8.0));
	clear_latvec( F_OFFSET(xxx2), EVENANDODD );
	restore_fermion_links_from_site(fn_links, MILC_PRECISION);
	fn = get_fm_links(fn_links);
     	grsource_imp( F_OFFSET(phi2), mass2, EVEN, fn[n-1]);

	/* update U's to middle of interval */
     	update_u(epsilon*nflavors2/8.0);
#endif

        /* do conjugate gradient to get (Madj M)inverse * phi */
	restore_fermion_links_from_site(fn_links, MILC_PRECISION);
	fn = get_fm_links(fn_links);
	if(n == 2){
	  iters += ks_congrad( F_OFFSET(phi1), F_OFFSET(xxx1), mass1,
		       niter, nrestart, rsqmin, MILC_PRECISION, EVEN, &final_rsq, fn[0] );
	  iters += ks_congrad( F_OFFSET(phi2), F_OFFSET(xxx2), mass2,
		       niter, nrestart, rsqmin, MILC_PRECISION, EVEN, &final_rsq, fn[1] );
	} else {
	iters += ks_congrad_two_src( F_OFFSET(phi1), F_OFFSET(phi2),
				     F_OFFSET(xxx1), F_OFFSET(xxx2),
				     mass1, mass2, niter, nrestart, rsqmin, 
				     MILC_PRECISION, EVEN, &final_rsq,
				     fn[0]);
	}

	dslash_site( F_OFFSET(xxx1), F_OFFSET(xxx1), ODD, fn[0]);
	dslash_site( F_OFFSET(xxx2), F_OFFSET(xxx2), ODD, fn[n-1]);
	/* now update H by full time interval */
    	update_h(epsilon);

#if 0
#ifdef HAVE_QIO
	{
	  char *filexml;
	  char recxml[] = "<?xml version=\"1.0\" encoding=\"UTF-8\"?><title>Test fermion force field</title>";
	  char ansfile[128];
	  char rootname[] = "fermion_force_dump";

	  /* Do this at the specified interval */
	  if(step%3 == 1){
	    
	    /* Construct a file name */
	    sprintf(ansfile,"%s%02d",rootname,step);
	    /* Dump the computed fermion force from the site structure */
	    filexml = create_QCDML();
	    save_color_matrix_scidac_from_site(ansfile, filexml, 
		       recxml, QIO_PARTFILE,  F_OFFSET(mom[0]), 4, MILC_PRECISION);
	    free_QCDML(filexml);
	  }
	}
#endif
#endif

    	/* update U's by half time step to get to even time */
    	update_u(epsilon*0.5);

        /* reunitarize the gauge field */
        reunitarize_ks();

    }	/* end loop over microcanonical steps */

#ifdef HMC_ALGORITHM
    /* find action */
    /* do conjugate gradient to get (Madj M)inverse * phi */
    restore_fermion_links_from_site(fn_links, MILC_PRECISION);
    fn = get_fm_links(fn_links);
    iters += ks_congrad( F_OFFSET(phi1), F_OFFSET(xxx1), mass1,
			 niter, nrestart, rsqmin, MILC_PRECISION, EVEN, 
			 &final_rsq, fn[0]);
    iters += ks_congrad( F_OFFSET(phi2), F_OFFSET(xxx2), mass2,
			 niter, nrestart, rsqmin, MILC_PRECISION, EVEN, 
			 &final_rsq, fn[n-1]);
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


/**********************************************************************/
/*   Accessor for string describing the option                        */
/**********************************************************************/
const char *ks_int_alg_opt_chr( void )
{
  return "INT_ALG_NEEDS_TO_BE_FIXED";
}

