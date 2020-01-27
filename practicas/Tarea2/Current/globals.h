#ifndef __GLOBALS_H
#define __GLOBALS_H

//#include "mpi.h"

#include "first_macro_1d.h"  /* Where global parameters are defined */
#include "structs_1d.h"      /* Where structures are defined */

struct globals {

    /* ---------------------------  rkc -------------------------------------*/

#ifdef IMEX

    struct field_array dv_tempF2;
    struct field_array dv_tempF3;
    struct field_array dv_tempF4;
    struct field_array dv_tempS1;
    struct field_array dv_tempS2;
    struct field_array dv_tempS3;
    struct field_array v_temp0;
    struct field_array v_temp1;

#else
    
    struct field_array dv_temp;        /* derivative temporal value         */
    struct field_array dv_sum;         /* sum values of deriv               */
    struct field_array v_temp;         /* intermediate value of v           */

#endif

    /* ---------------------------- FF --------------------------------------*/
    struct deriv_array dfields;
    struct aux_array auxfields;

};

extern struct globals globals;

#endif
