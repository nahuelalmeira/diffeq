#ifndef _DVDEF
#define _DVDEF

#ifdef HAVE_IEEEFP_H
#include <ieeefp.h>
#endif

#include <stdio.h>
#include "v_types.h"

extern     double       dv_dmin(double d1,double d2);
extern     double       dv_dmax(double d1,double d2);

extern     DVEC         make_DVEC(r_int n);
extern     void         free_DVEC(r_DVEC v);
extern     DVEC         make_DVEC_p(r_int n);
extern     void         free_DVEC_p(r_DVEC v);

extern     DVECN        make_DVECN(int n);
extern     void         free_DVECN(DVECN the);

extern     DIVEC        make_DIVEC(r_int n);
extern     void         free_DIVEC(DIVEC v);

extern     double       dvmin(r_DVEC v,r_int n);
extern     double       dvmax(r_DVEC v,r_int n);
extern     double       dvabsmin(r_DVEC v,r_int n);
extern     double       dvabsmax(r_DVEC v,r_int n);
extern     double       dxymin_range(r_DVEC x,r_DVEC y,r_int n,
                           r_double xmin,r_double xmax);
extern     double       dxymax_range(r_DVEC x,r_DVEC y,r_int n,
                           r_double xmin,r_double xmax);
extern     void         dvprint(DVEC v,int n,char *s);
extern     void         divprint(DIVEC v,int n,char *s);
extern     void         dvprint_inc(DVEC v,int n,int inc,char *s);
extern     void         dvprint_range(DVEC v,int n,char *s);
extern     DVEC         Dviota(r_int n);
extern     DVEC         DViota(r_DVEC v,r_int n);
extern     DVEC         Dvcopy(r_DVEC v1,r_int n);
extern     DVEC         DVcopy(r_DVEC v1,r_DVEC v2,r_int n);

extern     void         dvva(r_DVEC v1,r_DVEC v2,r_DVEC v3,r_int n);
extern     void         dvvm(r_DVEC v1,r_DVEC v2,r_DVEC v3,r_int n);
extern     void         dvvs(r_DVEC v1,r_DVEC v2,r_DVEC v3,r_int n);
extern     void         dvsm(r_DVEC v1,r_double s1,r_DVEC v2,int n);
extern     void         dvsa(r_DVEC v1,r_double s1,r_DVEC v2,int n);

extern     DVEC         DVva(r_DVEC v1,r_DVEC v2,r_DVEC v3,r_int n);
extern     DVEC         DVvm(r_DVEC v1,r_DVEC v2,r_DVEC v3,r_int n);
extern     DVEC         DVvs(r_DVEC v1,r_DVEC v2,r_DVEC v3,r_int n);
extern     DVEC         DVsm(r_DVEC v1,r_double s1,r_DVEC v2,int n);
extern     DVEC         DVsa(r_DVEC v1,r_double s1,r_DVEC v2,int n);

extern     DVEC         DVrand(r_DVEC v,r_double v_min,r_double v_max,r_int n);

extern     int          DVscanf(r_DVEC v,r_int n);
extern     int          DVprintf(r_DVEC v,r_int n);
extern     int          DVVscanf(r_DVEC v1,r_DVEC v2,r_int n);
extern     int          DVVprintf(r_DVEC v1,r_DVEC v2,r_int n);

extern     void         dvls(r_DVEC v1,double s1,int n);
extern     void         dvav(r_DVEC v1,r_DVEC v2,int n);
extern     int          idvnz(DVEC v,int n);
extern     void         dvramp(DVEC v,double v0,double dv,int n);
extern     DVEC         Dvramp(double v0,double dv,int n);
extern     void         dvumsh(DVEC v,int n,double v0,double vnm1);

extern     double       dvsum(DVEC v,int n);
extern     double       dvvsum(DVEC v1,DVEC v2,int n);

extern     void         dvdd01(DVEC v1,DVEC v2,double h,int n);

extern     void         dvdump(DVEC v,int n,char *s);
extern     void         dvpdmp(DVEC v1,DVEC v2,int n,char *s);
extern     void         dvfdump(FILE *fp,DVEC v,int n,char *s);

extern     double       drange(double xmin,double x,double xmax,double fuzz);
extern     void         generate_code(char *f_name,char *r_name,char *sub_name,
                                      DVEC c,int nc,double xmin,double xmax);
extern     void         generate_c_data(FILE *code_file,DVEC c,int nc,int nseg);
extern     void         generate_xseg_data(FILE *code_file,DVEC xseg,int nseg);
extern     void         generate_xseg_data_1(FILE *code_file,DVEC xseg,
                           int fi_st,int fi_fin);
extern     void         old_generate_xseg_data(FILE *code_file,
                            DVEC xseg,int nseg);
extern     void         generate_multi_segment_code(char *f_name,char *r_name,
                           char *sub_name,DVEC c,int nc,DVEC xseg,int  nseg);
extern     void         generate_multi_segment_code_x(char *f_name,
                           char *r_name,char *sub_name,DVEC c,int nc,DVEC xseg,
                           int nseg);
extern     void         generate_multi_segment_code_dn(char *f_name, 
                           char *r_name, char *sub_name, 
                           DVEC c, int nc,DVEC xseg, int nseg);
extern     void         generate_dc_and_c_data(FILE *code_file, DVEC c, int nc, 
                           int nseg);

extern     void         dvfout(FILE *fout,DVEC v,int n);

extern     double       dvmindv(DVEC v,int n);

extern     void         dvsortup(DVEC v, int n);
extern     void         dvsortdown(DVEC v, int n);
extern     int          sort_up(double *a, double *b); 
extern     int          sort_down(double *a,double *b);
extern     void         dvsortupuniq_(DVEC vin,int *pnin,DVEC vout,
                           int *pnout,double *pfuzz);
extern     void         dvsortupuniq(DVEC vin, int nin, DVEC vout, 
                           int *pnout, double fuzz);
extern     void         dvsortdownuniq(DVEC vin, int nin, DVEC vout, 
                           int *pnout, double fuzz);
extern     int          dvlookup_linear(DVEC v,int n,double vkey,double fuzz);
extern     int          eq_fuzz(double x1,double x2,double fuzz);

extern     int          dvlsindex(DVEC v,int n,double key,double fuzz);
extern     IVEC         make_IVEC_dv(r_int n);
extern     IVEC         Dvlsindex(DVEC v,int n,DVEC key,int nkey,double fuzz);

extern     void         dvvmergeup(DVEC v1,int n1,DVEC v2,int n2,
                           DVEC v3,int *pn3,double fuzz);
extern     void         dvvmergedown(DVEC v1,int n1,DVEC v2,int n2,
                           DVEC v3,int *pn3,double fuzz);

extern     void         dvmrg(DVEC v1,DVEC v2,DVEC v3,int n);
extern     void         dvumrg(DVEC v1,DVEC v2,DVEC v3,int n);

#define    SORT_UP      0
#define    SORT_DOWN    1

extern     int          dv2sort(DVEC v1,DVEC v2,int n,int code);

extern     int          sort_up_d2(PD2 el1,PD2 el2);
extern     int          sort_down_d2(PD2 el1,PD2 el2);

extern     DVEC         dvdel_el(DVEC v,int n,int i);
extern     int          dv2insert(DVEC *pv1,DVEC *pv2,int n,
                           double new1,double new2,int code);
extern     int          dv2insert_nosort(DVEC *pv1,DVEC *pv2,int n,
                           double new1,double new2,int i);

extern     int          ixdvmin(DVEC v,int n);
extern     int          ixdvmax(DVEC v,int n);
extern     int          ixdvabsmin(DVEC v,int n);
extern     int          ixdvabsmax(DVEC v,int n);

extern     int          ixdvnearest(DVEC v,int n,double key);
extern     int          ixdv2nearest(DVEC v1,DVEC v2,int n,
                           double key1,double key2,double asp1by2);

extern     double       dvhash(DVEC v,int n);

extern     void         dvfapl(DVEC v1,DVEC v2,PFD f,int n);
extern     DVEC         Dvfapl(DVEC v1,PFD f,int n);

extern     DVEC         dvfscanf(FILE *fp,DVEC v,int *pn);
extern     void         dvfprintf(FILE *fp,DVEC v,int n);

extern     int          dvget(char *fname,DVEC v,int *pn);
extern     DVEC         Dvget(char *fname,int *pn);
extern     int          dvput(char *fname,DVEC v,int n);

extern     void         dvvfscanf(FILE *fp,DVEC v1,DVEC v2,int *pn);
extern     int          dvvget(char *fname,DVEC v1,DVEC v2,int *pn);
extern     int          dvvput(char *fname,DVEC v1,DVEC v2,int n);  
extern     void         dvvfprintf(FILE *fp,DVEC v1,DVEC v2,int n);
extern     DVEC         Dvfscanf(FILE *fp,int *pn);
extern     PDVEC        Dvvfscanf(FILE *fp,int *pn);
extern     PDVEC        Dvvget(char *fname,int *pn);

extern     void         dvvvfscanf(FILE *fp,DVEC v1,DVEC v2,DVEC v3,int *pn);

extern     PDVECN       PDVECNget(char *fname);
extern     PDVECN       free_and_null_PDVECN(PDVECN p);

extern     void         PDVECNfdump(FILE *fp,PDVECN the,char *s);

extern     void         dvpyth(DVEC v1,DVEC v2,DVEC v3,int n);
extern     void         dvavvm(DVEC v1,DVEC v2,DVEC v3,DVEC v4,int n);

extern     DVEC         dvappend(DVEC v,int n,double i);
extern     DVECN        dvecn_append(DVECN the,double i);
extern     int          dv_between(double val1,double val,double val2);
extern     DVEC         Dv2perp(DVEC vx,DVEC vy,int n,
                           double x,double y,double aspect);
extern     int          ixdv2perp(DVEC vx,DVEC vy,int n,
                           double x,double y,double aspect);

extern     SVEC         DVSVcopy(DVEC vin,SVEC vout,int n);

extern     PDVECN       PDVECN_interp(PDVECN from,PDVECN to,
                                      double vs,double vf,int intord);
extern     PDVECN       PDVECN_lintseg(PDVECN in,double xscale,double yscale,
                                       double mindseg,double maxdseg);

extern     IVEC         ivmakemaskdv(DVEC vin,int nin,PFI_D logic);
extern     DVEC         dvreduce(DVEC vin,IVEC mask,int nin,int *pnout);
extern     int          isdvramp(DVEC v,int n,double fuzz);

extern     int          limin(int i1,int i2);
extern     int          limax(int i1,int i2);

extern     void         dvinqn(DVEC v,DVEC x,DVEC vbar,DVEC xbar,int n,int nbar,
                               double vs, double vf,int nintrp);
extern     double       flipn(double xbar,DVEC x,DVEC y,int n);

extern     double       fdvx0(DVEC x,DVEC y,int n,int jst,int *pjx);
extern     void         dvx0(DVEC x,DVEC y,int n,DVEC x0,int *pn0);

extern     int          isdvfinite(DVEC v,int n);

#define     MDVFAIL(v,n)    !(v = make_DVEC(n))

#include  "dveclib.h"
#include  "dveclib1.h"

#endif
