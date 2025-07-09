#ifndef _PARAMS_H
#define _PARAMS_H

#include "../include/macros.h"  /* For MAXFILENAME */
#include "defines.h"
#include "../include/imp_ferm_links.h"

/* structure for passing simulation parameters to each node */
typedef struct {
	int stopflag;   /* 1 if it is time to stop */
    /* INITIALIZATION PARAMETERS */
	int nx,ny,nz,nt;  /* lattice dimensions */
#ifdef FIX_NODE_GEOM
  int node_geometry[4];  /* Specifies fixed "nsquares" (i.e. 4D
			    hypercubes) for the compute nodes in each
			    coordinate direction.  Must be divisors of
			    the lattice dimension */
#ifdef FIX_IONODE_GEOM
  int ionode_geometry[4]; /* Specifies fixed "nsquares" for I/O
			     partitions in each coordinate direction,
			     one I/O node for each square.  The I/O
			     node is at the origin of the square.
			     Must be divisors of the node_geometry. */
#endif
#endif
	int iseed;	/* for random numbers */
#ifdef ONEMASS
	int nflavors;	/* the number of flavors */
#else
	int nflavors1;	/* the number of flavors of first type */
	int nflavors2;	/* the number of flavors of second type */
#endif
    /*  REPEATING BLOCK */
	int warms;	/* the number of warmup trajectories */
	int trajecs;	/* the number of real trajectories */
	int steps;	/* number of steps for updating */
	int propinterval;     /* number of trajectories between measurements */
	Real beta;      /* gauge coupling */
#ifdef ONEMASS
	Real mass;      /*  quark mass */
	#ifdef HASENBUSCH
		Real hmass;
	#endif
#else
	Real mass1,mass2; /*  quark masses */
        Real naik_term_epsilon2;   /* Naik term parameter for 2nd mass */
#endif
	Real u0; /* tadpole parameter */
	int aniter; 	/* maximum number of c.g. iterations */
	int fniter; 	/* maximum number of c.g. iterations */
    int anrestart;   /* maximum number of c.g. restarts */
	int fnrestart;   /* maximum number of c.g. restarts */
        ks_eigen_param eigen_param; /* Parameters for eigensolver. Not used for HMC */
        int npbp_reps_in;   /* Number of random sources */
        int prec_pbp;       /* Precision of pbp measurements */
	Real rsqprop;  /* for deciding on convergence */
	Real arsqmin,frsqmin;  /* for deciding on convergence */
	Real epsilon;	/* time step */
        char spectrum_request[MAX_SPECTRUM_REQUEST];   /* request list for spectral measurements */
        int source_start, source_inc, n_sources; /* source time and increment */
        int spectrum_multimom_nmasses; 
        Real spectrum_multimom_low_mass;
        Real spectrum_multimom_mass_step;
        int fpi_nmasses;
        Real fpi_mass[MAX_FPI_NMASSES];
	int startflag;  /* what to do for beginning lattice */
	int saveflag;   /* what to do with lattice at end */
	char startfile[MAXFILENAME],savefile[MAXFILENAME];
	char stringLFN[MAXFILENAME];  /** ILDG LFN if applicable ***/
}  params;

#endif /* _PARAMS_H */
