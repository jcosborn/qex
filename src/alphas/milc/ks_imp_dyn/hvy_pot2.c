/****** hvy_pot2.c  -- ******************/
/* MIMD version 7 */
/* Heavy quark potential
* MIMD version 7
/* NOT MAINTAINED.  TEST BEFORE USE! */
/* DT 6/10/97
* DT 7/17/97 modified to do all displacements
* Evaluate in different spatial directions, to check rotational
* invariance.  Gauge fix to Coulomb gauge to do spatial segments.
*
* Argument tells whether to use ordinary links or fat links.
* Gauge must be fixed to Coulomb gauge first.
*/

#include "ks_imp_includes.h"	/* definitions files and prototypes */
void printpath( int *path, int length );
void shiftmat( field_offset src, field_offset dest, int dir );

void hvy_pot( field_offset links, int max_t, int max_x ) {
    register int i,j,k;
    register site *s;
    int t_dist, x_dist, y_dist, z_dist;
    double wloop;
    msg_tag *mtag0;
    field_offset oldmat, newmat, tt;

    node0_printf("hvy_pot(): MAX_T = %d, MAX_X = %d\n",max_t,max_x);

    rephase( OFF );

    /* Use tempmat1 to construct t-direction path from each point */
    for( t_dist=1; t_dist<= max_t; t_dist ++){
	if(t_dist==1 ){
	    FORALLSITES(i,s){
		su3mat_copy( &(((su3_matrix *)(F_PT(s,links)))[TUP]),
		    &(s->tempmat1) );
	    }
	}
	else{
	    mtag0 = start_gather_site( F_OFFSET(tempmat1), sizeof(su3_matrix),
	        TUP, EVENANDODD, gen_pt[0] );
	    wait_gather(mtag0);
	    FORALLSITES(i,s){
	 	mult_su3_nn( &(((su3_matrix *)(F_PT(s,links)))[TUP]),
		    (su3_matrix *)gen_pt[0][i], &(s->staple) );
	    }
	    cleanup_gather( mtag0 );
	    FORALLSITES(i,s){
	 	su3mat_copy( &(s->staple), &(s->tempmat1) );
	    }
	}
	/* now tempmat1 is path of length t_dist in TUP direction */
	oldmat= F_OFFSET(tempmat2);
	newmat= F_OFFSET(staple);	/* will switch these two */

	for( x_dist=0; x_dist<=max_x; x_dist++ ){
	    /**for( y_dist=0; y_dist<=x_dist; y_dist+= x_dist>0?x_dist:1 ){**/
	    for( y_dist=0; y_dist<=max_x; y_dist+= 1 ){
		/* now gather from spatial dirs, compute products of paths */
		FORALLSITES(i,s){
		    su3mat_copy( &(s->tempmat1), (su3_matrix *)F_PT(s,oldmat) );
		}
		for(i=0;i<x_dist;i++){
		    shiftmat( oldmat, newmat, XUP );
		    tt=oldmat; oldmat=newmat; newmat=tt;
		}
		for(i=0;i<y_dist;i++){
		    shiftmat( oldmat, newmat, YUP );
		    tt=oldmat; oldmat=newmat; newmat=tt;
		}
		for( z_dist=0; z_dist<=max_x; z_dist++ ){
		    /* evaluate potential at this separation */
		    wloop = 0.0;
		    FORALLSITES(i,s){
			wloop += (double)realtrace_su3( &(s->tempmat1),
			    (su3_matrix *)F_PT(s,oldmat) );
		    }
		    g_doublesum( &wloop );
		    node0_printf("POT_LOOP: %d %d %d %d \t%e\n",
			x_dist, y_dist, z_dist, t_dist, wloop/volume );

		    /* as we increment z, shift in z direction */
		    shiftmat( oldmat, newmat, ZUP );
		    tt=oldmat; oldmat=newmat; newmat=tt;
		} /* z distance */
	    } /* y distance */
	} /*x distance */
    } /*t_dist*/
    rephase( ON );
} /* hvy_pot() */

/* shift, without parallel transport, a matrix from direction "dir" */
void shiftmat( field_offset src, field_offset dest, int dir ){
    register int i;
    register site *s;
    msg_tag *mtag;
    mtag = start_gather_site( src, sizeof(su3_matrix),
        dir, EVENANDODD, gen_pt[0] );
    wait_gather(mtag);
    FORALLSITES(i,s){
        su3mat_copy( (su3_matrix *)gen_pt[0][i], (su3_matrix *)F_PT(s,dest) );
    }
    cleanup_gather( mtag );
}
