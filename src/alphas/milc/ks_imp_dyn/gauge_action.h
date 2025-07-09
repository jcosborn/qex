/* Symanzik 1-loop gauge action: 1x1 + 1x2 + 1x1x1 with quark loops
   included */
#ifdef GAUGE_ACTION_PART1
#define SYMANZIK_1LOOP
/* defines NREPS NLOOP MAX_LENGTH MAX_NUM */
#undef MAX_LENGTH
#undef MAX_NUM
#define NREPS 1
#define NLOOP 3
#define MAX_LENGTH 6
#define MAX_NUM 16
#endif

#ifdef GAUGE_ACTION_PART2
    static int loop_ind[NLOOP][MAX_LENGTH] = {
    { XUP, YUP, XDOWN, YDOWN, NODIR, NODIR },
    { XUP, XUP, YUP, XDOWN, XDOWN , YDOWN},
    { XUP, YUP, ZUP, XDOWN, YDOWN , ZDOWN},
    };
    static int loop_length_in[NLOOP] = {4,6,6};

    for(j=0;j<NLOOP;j++){
	loop_num[j] = 0;
	loop_length[j] = loop_length_in[j];
	for(i=0;i<NREPS;i++){
	    loop_coeff[j][i] = 0.0;
	}
    }

    /* Loop coefficients from Urs */
    /* For the Asqtad flavor dependence, see Zh. Hao, G.M von Hippel, RR
       Horgan, QJ Mason, and HD Trottier, hep-lat/0705.4660v2 */
    // HISQ coefficients: Ron Horgan private communication
    loop_coeff[0][0]= 1.0;
    loop_coeff[1][0]=  -1.00/(20.0*u0*u0) * 
      ( 1.00 - (0.6264 - 1.1746*total_dyn_flavors)*log(u0) );
    loop_coeff[2][0]=  1.00/(u0*u0) * 
      (0.0433 - 0.0156*total_dyn_flavors) * log(u0); 
    strcpy(gauge_action_description,"\"Symanzik 1x1 + 1x2 + 1x1x1 action with HISQ quark loops\"");
    node0_printf("Symanzik 1x1 + 1x2 + 1x1x1 action with HISQ quark loops\n");

#endif
