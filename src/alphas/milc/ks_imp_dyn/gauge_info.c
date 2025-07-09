/*********************** gauge_info.c *************************/
/* MIMD version 7 */

/* For ks_imp_dyn */

/* Application-dependent routine for writing gauge info file */
/* This file is an ASCII companion to the gauge configuration file
   and contains information about the action used to generate it.
   This information is consistently written in the pattern

       keyword  value

   or

       keyword[n] value1 value2 ... valuen

   where n is an integer.

   To maintain a semblance of consistency, the possible keywords are
   listed in io_lat.h.  Add more as the need arises, but be sure
   to notify the rest of the collaboration.

   */

#include "ks_imp_includes.h"
#include <quark_action.h>
#ifdef HAVE_QIO
#include <qio.h>
#endif

/*---------------------------------------------------------------------------*/
/* This routine writes the ASCII info file.  It is called from one of
   the lattice output routines in io_lat4.c.*/

extern char gauge_action_description[128]; /* in gauge_stuff.c */
extern int gauge_action_nloops,gauge_action_nreps;
void write_appl_gauge_info(FILE *fp, gauge_file *gf)
{
  Real myssplaq = g_ssplaq;  /* Precision conversion */
  Real mystplaq = g_stplaq;  /* Precision conversion */
  Real nersc_linktr = linktrsum.real/3.;  /* Convention and precision */

  /* Write generic information */
  write_generic_gauge_info(fp, gf);

  /* The rest are optional */

  write_gauge_info_item(fp,"action.description","%s",
			"\"Gauge plus fermion\"",0,0);

  write_gauge_info_item(fp,"gauge.description","%s",
			gauge_action_description,0,0);
  write_gauge_info_item(fp,"gauge.nloops","%d",(char *)&gauge_action_nloops,0,0);
  write_gauge_info_item(fp,"gauge.nreps","%d",(char *)&gauge_action_nreps,0,0);
  write_gauge_info_item(fp,"gauge.beta11","%f",(char *)&beta,0,0);
  write_gauge_info_item(fp,"gauge.tadpole.u0","%f",(char *)&u0,0,0);
  write_gauge_info_item(fp,"gauge.ssplaq","%f",(char *)&myssplaq,0,0);
  write_gauge_info_item(fp,"gauge.stplaq","%f",(char *)&mystplaq,0,0);
  write_gauge_info_item(fp,"gauge.nersc_linktr","%e",
			(char *)&(nersc_linktr),0,0);
  write_gauge_info_item(fp,"gauge.nersc_checksum","%u",
			(char *)&(nersc_checksum),0,0);
  write_gauge_info_item(fp,"quark.description","%s",QUARK_ACTION_DESCRIPTION,0,0);
#ifdef ONEMASS
  write_gauge_info_item(fp,"quark.flavors","%d",(char *)&nflavors,0,0);
  write_gauge_info_item(fp,"quark.mass","%f",(char *)&mass,0,0);
#else
  write_gauge_info_item(fp,"quark.flavors1","%d",(char *)&nflavors1,0,0);
  write_gauge_info_item(fp,"quark.flavors2","%d",(char *)&nflavors2,0,0);
  write_gauge_info_item(fp,"quark.mass1","%f",(char *)&mass1,0,0);
  write_gauge_info_item(fp,"quark.mass2","%f",(char *)&mass2,0,0);
#endif
}

/* Print string formatting with digits based on intended precision */
void print_prec(char string[], size_t n, Real value, int prec){
  if(prec == 1){
    snprintf(string,n,"%.6e",value);  /* single precision */
  }
  else
    snprintf(string,n,"%.15e",value); /* double precision */
}

#define INFOSTRING_MAX 4096
char *create_MILC_info(){

  size_t bytes = 0;
  char *info = (char *)malloc(INFOSTRING_MAX);
  size_t max = INFOSTRING_MAX;
  Real myssplaq = g_ssplaq;  /* Precision conversion */
  Real mystplaq = g_stplaq;  /* Precision conversion */
  Real nersc_linktr = linktrsum.real/3.;  /* Convention and precision */

  sprint_gauge_info_item(info+bytes, max-bytes,"action.description","%s",
			"\"Gauge plus fermion\"",0,0);
  bytes = strlen(info);
  sprint_gauge_info_item(info+bytes, max-bytes,"gauge.description","%s",
			gauge_action_description,0,0);
  bytes = strlen(info);
  sprint_gauge_info_item(info+bytes, max-bytes,"gauge.nloops","%d",
			 (char *)&gauge_action_nloops,0,0);
  bytes = strlen(info);
  sprint_gauge_info_item(info+bytes, max-bytes,"gauge.nreps","%d",
			 (char *)&gauge_action_nreps,0,0);
  bytes = strlen(info);
  sprint_gauge_info_item(info+bytes, max-bytes,"gauge.beta11","%f",
			 (char *)&beta,0,0);
  bytes = strlen(info);
  sprint_gauge_info_item(info+bytes, max-bytes,"gauge.tadpole.u0","%f",
			 (char *)&u0,0,0);
  bytes = strlen(info);
  sprint_gauge_info_item(info+bytes, max-bytes,"gauge.ssplaq","%f",
			 (char *)&myssplaq,0,0);
  bytes = strlen(info);
  sprint_gauge_info_item(info+bytes, max-bytes,"gauge.stplaq","%f",
			 (char *)&mystplaq,0,0);
  bytes = strlen(info);
  sprint_gauge_info_item(info+bytes, max-bytes,"gauge.nersc_linktr","%e",
			 (char *)&nersc_linktr,0,0);
  bytes = strlen(info);
  sprint_gauge_info_item(info+bytes, max-bytes,"gauge.nersc_checksum","%u",
			 (char *)&nersc_checksum,0,0);
  
  bytes = strlen(info);
  sprint_gauge_info_item(info+bytes, max-bytes,"quark.description","%s",
			 QUARK_ACTION_DESCRIPTION,0,0);
  bytes = strlen(info);
#ifdef ONEMASS
  sprint_gauge_info_item(info+bytes, max-bytes,"quark.flavors","%d",
			 (char *)&nflavors,0,0);
  bytes = strlen(info);
  sprint_gauge_info_item(info+bytes, max-bytes,"quark.mass","%f",
			 (char *)&mass,0,0);
  bytes = strlen(info);
#else
  sprint_gauge_info_item(info+bytes, max-bytes,"quark.flavors1","%d",
			 (char *)&nflavors1,0,0);
  bytes = strlen(info);
  sprint_gauge_info_item(info+bytes, max-bytes,"quark.flavors2","%d",
			 (char *)&nflavors2,0,0);
  bytes = strlen(info);
  sprint_gauge_info_item(info+bytes, max-bytes,"quark.mass1","%f",
			 (char *)&mass1,0,0);
  bytes = strlen(info);
  sprint_gauge_info_item(info+bytes, max-bytes,"quark.mass2","%f",
			 (char *)&mass2,0,0);
  bytes = strlen(info);
#endif  

  return info;
}

void destroy_MILC_info(char *info){
  if(info != NULL)
    free(info);
}


#ifdef HAVE_QIO
static QIO_String *xml_record;

/* Follow USQCD style for record XML */

char *create_QCDML(){
  QIO_USQCDLatticeInfo *record_info;
  char *milc_info;
  char plaqstring[32];
  char linktrstring[32];
  
  /* Build the components of the USQCD info string */
  print_prec(plaqstring, 32, (g_ssplaq+g_stplaq)/6., 1);
  print_prec(linktrstring, 32, linktrsum.real/3., 1);
  milc_info = create_MILC_info();  /* The MILC info data as a string */

  /* Stuff the data structure */
  record_info = QIO_create_usqcd_lattice_info(plaqstring, linktrstring, 
					      milc_info);

  destroy_MILC_info(milc_info);

  /* Convert the data structure to a QIO character string */
  xml_record = QIO_string_create();
  QIO_encode_usqcd_lattice_info(xml_record, record_info);
  QIO_destroy_usqcd_lattice_info(record_info);

  /* Return a pointer to the actual character string */
  return QIO_string_ptr(xml_record);
}

void free_QCDML(char *info){
  QIO_string_destroy(xml_record);
}

#endif
