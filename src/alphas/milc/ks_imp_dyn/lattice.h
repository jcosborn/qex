#ifndef _LATTICE_H
#define _LATTICE_H
/****************************** lattice.h ********************************/

/* include file for MIMD version 7
   This file defines global scalars and the fields in the lattice.

   Directory for dynamical improved KS action.  Allow:
	arbitrary paths in quark action 
	arbitrary paths in gauge action (eg Symanzik imp.)

   If "FN" is defined,
     Includes storage for Naik improvement (longlink[4], templongvec[4],
     gen_pt[16], etc.
*/

#include "defines.h"
#include "../include/generic_quark_types.h"
#include "../include/macros.h"    /* For MAXFILENAME */
#include "../include/io_lat.h"    /* For gauge_file */
#include "params.h"
#include "../include/generic_ks.h" /* For fn_links and ks_act_paths */
#include "../include/fermion_links.h"

/* Begin definition of site structure */

#include "../include/su3.h"
#include "../include/random.h"   /* For double_prn */

/* The lattice is an array of sites.  */
#define MOM_SITE   /* If there is a mom member of the site struct */
typedef struct {
    /* The first part is standard to all programs */
	/* coordinates of this site */
	short x,y,z,t;
	/* is it even or odd? */
	char parity;
	/* my index in the array */
	uint32_t index;
#ifdef SITERAND
	/* The state information for a random number generator */
	double_prn site_prn;
	/* align to double word boundary (kludge for Intel compiler) */
	int space1;
#endif

/* ------------------------------------------------------------ */
/*   Now come the physical fields, program dependent            */
/* ------------------------------------------------------------ */
	/* gauge field */
	su3_matrix link[4] ALIGNMENT;	/* the fundamental field */
#if  defined(HYBRIDS)
	su3_matrix tmplink[4];	/* three link straight paths */
#endif

#ifdef HMC_ALGORITHM
 	su3_matrix old_link[4];
	/* For accept/reject */
#endif

	/* antihermitian momentum matrices in each direction */
 	anti_hermitmat mom[4] ALIGNMENT;

	/* The Kogut-Susskind phases, which have been absorbed into 
		the matrices.  Also the antiperiodic boundary conditions.  */
 	Real phase[4];

	/* 3 element complex vectors */
#ifdef ONEMASS
 	su3_vector phi;		/* Gaussian random source vector */
 	su3_vector xxx;		/* solution vector = Kinverse * phi */
	su3_vector psi;
	#ifdef HASENBUSCH
	  su3_vector hphi;
	  su3_vector hxxx;
	  su3_vector hxxxr;
	  su3_vector hpsi;
	#endif
#ifdef HMC_ALGORITHM
  su3_vector old_xxx;	/* For predicting next xxx */
#endif
#else
 	su3_vector phi1;	/* Gaussian random source vector, mass1 */
 	su3_vector phi2;	/* Gaussian random source vector, mass2 */
 	su3_vector xxx1;	/* solution vector = Kinverse * phi, mass1 */
 	su3_vector xxx2;	/* solution vector = Kinverse * phi, mass2 */
#endif
 	su3_vector resid;	/* conjugate gradient residual vector */
 	su3_vector cg_p;	/* conjugate gradient change vector */
 	su3_vector ttt;		/* temporary vector, for K*ppp */
 	su3_vector g_rand;	/* Gaussian random vector*/
	/* Use trick of combining xxx=D^adj D)^(-1) on even sites with
	   Dslash times this on odd sites when computing fermion force */
	
#ifdef HYBRIDS
        su3_matrix field_strength[6];
#endif
#ifdef SPECTRUM
        su3_matrix tempmat1,staple;
	su3_vector propmat[3];	/* For three source colors */
	su3_vector propmat2[3];	/* nl_spectrum() */
	su3_matrix tempmat2;
	/* for spectrum_imp() */
	/**su3_vector quark_source, quark_prop, anti_prop;**/
#define quark_source propmat2[0]
#define quark_prop propmat2[1]
#define anti_prop propmat2[2]
#endif

	/* temporary vectors and matrices */
	su3_vector tempvec[4];	/* One for each direction */
#ifdef FN
	su3_vector templongvec[4];	/* One for each direction */
        su3_vector templongv1;
#endif
#ifdef DM_DU0
	su3_vector dMdu_x;	/* temp vector for <psi-bar(dM/du0)psi>
				   calculation in 'f_meas.c' */
#endif
#ifdef NPBP_REPS
 	su3_vector M_inv;	/* temp vector for M^{-1} g_rand */
#endif
#if defined(CHEM_POT) || defined(D_CHEM_POT)
 	su3_vector dM_M_inv;	/* temp vector for dM/dmu M^{-1} g_rand */
        su3_vector deriv[6];
#endif
} site;

/* End definition of site structure */

/* Definition of globals */

#ifdef CONTROL
#define EXTERN 
#else
#define EXTERN extern
#endif

/* The following are global scalars 
   beta is overall gauge coupling factor
   mass, mass1 and mass2 are quark masses
   u0 is tadpole improvement factor, perhaps (plaq/3)^(1/4)
*/
EXTERN	int nx,ny,nz,nt;	/* lattice dimensions */
EXTERN  size_t volume;		/* volume of lattice = nx*ny*nz*nt */
#ifdef FIX_NODE_GEOM
EXTERN  int node_geometry[4];  /* Specifies fixed "nsquares" (i.e. 4D
			    hypercubes) for the compute nodes in each
			    coordinate direction.  Must be divisors of
			    the lattice dimensions */
#ifdef FIX_IONODE_GEOM
EXTERN int ionode_geometry[4]; /* Specifies fixed "nsquares" for I/O
			     partitions in each coordinate direction,
			     one I/O node for each square.  The I/O
			     node is at the origin of the square.
			     Must be divisors of the node_geometry. */
#endif
#endif
EXTERN  params param;
EXTERN	uint32_t iseed;		/* random number seed */
EXTERN	int warms,trajecs,steps,niter,nrestart,propinterval;
EXTERN  int aniter,fniter,anrestart,fnrestart;
EXTERN  int npbp_reps_in;
EXTERN  int prec_pbp;  /* Precisiong of pbp measurements */
EXTERN  int dyn_flavors[MAX_DYN_MASSES]; 
#ifdef ONEMASS
EXTERN  int nflavors;
#else
EXTERN	int nflavors1,nflavors2;  /* number of flavors of types 1 and 2 */
#endif
EXTERN  int nlight_flavors;
EXTERN	Real epsilon;
EXTERN  Real beta,u0;
EXTERN  int n_dyn_masses; // number of dynamical masses
#ifdef ONEMASS
	EXTERN  Real mass;
	#ifdef HASENBUSCH
		EXTERN Real hmass;
	#endif
#else
EXTERN  Real mass1,mass2;
#endif
EXTERN  Real naik_term_epsilon2;
EXTERN	Real arsqmin,frsqmin;
EXTERN  Real rsqmin,rsqprop;
EXTERN	int startflag;	/* beginning lattice: CONTINUE, RELOAD, RELOAD_BINARY,
			   RELOAD_CHECKPOINT, FRESH */
EXTERN	int saveflag;	/* do with lattice: FORGET, SAVE, SAVE_BINARY,
			   SAVE_CHECKPOINT */
EXTERN	char startfile[MAXFILENAME],savefile[MAXFILENAME];
EXTERN  double g_ssplaq, g_stplaq;
EXTERN  double_complex linktrsum;
EXTERN  u_int32type nersc_checksum;
EXTERN  char stringLFN[MAXFILENAME];  /** ILDG LFN if applicable **/
EXTERN	int total_iters;
EXTERN	int hisq_svd_counter;
EXTERN	int hisq_force_filter_counter;
EXTERN	int hyphisq_svd_counter;
EXTERN  int hypisq_force_filter_counter;
EXTERN  int phases_in; /* 1 if KS and BC phases absorbed into matrices */
EXTERN  int source_start, source_inc, n_sources;
        /* source time, increment for it, and number of source slices */
EXTERN  char spectrum_request[MAX_SPECTRUM_REQUEST]; /* request list for spectral measurements */
/* parameters for spectrum_multimom */
EXTERN  int spectrum_multimom_nmasses;
EXTERN  Real spectrum_multimom_low_mass;
EXTERN  Real spectrum_multimom_mass_step;
/* parameters for fpi */
EXTERN  int fpi_nmasses;
EXTERN  Real fpi_mass[MAX_FPI_NMASSES];

/* Some of these global variables are node dependent */
/* They are set in "make_lattice()" */
EXTERN	size_t sites_on_node;		/* number of sites on this node */
EXTERN	size_t even_sites_on_node;	/* number of even sites on this node */
EXTERN	size_t odd_sites_on_node;	/* number of odd sites on this node */
EXTERN	int number_of_nodes;	/* number of nodes in use */
EXTERN  int this_node;		/* node number of this node */

EXTERN gauge_file *startlat_p;
EXTERN char hostname[128];

/* Each node maintains a structure with the pseudorandom number
   generator state */
EXTERN double_prn node_prn ;


/* The lattice is a single global variable - (actually this is the
   part of the lattice on this node) */
EXTERN Real boundary_phase[4];
EXTERN site *lattice;

/* Vectors for addressing */
/* Generic pointers, for gather routines */
#define N_POINTERS 16
EXTERN char ** gen_pt[N_POINTERS];

/* Storage for definition of the quark action */
EXTERN fermion_links_t        *fn_links;

/* For HISQ or HYPISQ operation */
EXTERN int n_order_naik_total;
EXTERN int n_pseudo_naik[MAX_N_PSEUDO];
EXTERN int n_orders_naik[MAX_N_PSEUDO];

/* For eigenpair calculation */
/* Not used for HMC */
EXTERN int Nvecs_tot;
EXTERN Real *eigVal; /* eigenvalues of D^dag D */
EXTERN su3_vector **eigVec; /* eigenvectors */
#endif /* _LATTICE_H */
