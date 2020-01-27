/*********************************************************************
* Here are the global definitions used by all the subroutines        *
*                                                                    *
*********************************************************************/

#ifndef FIRST_MACRO_1D_H
#define FIRST_MACRO_1D_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
/* #include <curses.h> */
/* #include <ncurses.h> */
#include <time.h>
#include <string.h>    /* functions to get data via web (stadin) */
/* #include <rfftw.h>  */    /*  needed for taking fft in derivs and main */

// #define DEBUG
// #define DEBUG_ADISCO

#define EXCISION
#undef EXCISION

#define DISSIPATION
#undef DISSIPATION

#undef STRICT


#define SIMPLE_BUMP
//#define SQUARE_BUMP
//#define BUMP_TEST

#define PERIODIC
//#define NO_LAST_POINT   // in order not to plot the last point when there are interfases


/* That is the number of variables that enter the system            */
#define N_FIELDS 2

/* The name of the fields */
#define U  0
#define V  1

#ifdef DISSIPATION
	//printf("dissipation not implemented!"); exit(0);
	#define N_DERIVS 4
	#define U_X 0
	#define V_X 1
	#define DISS_U 2
	#define DISS_V 3
#else
	#define N_DERIVS 2
	#define U_X 0
	#define V_X 1
#endif // DISSIPATION

#define N_AUX 4
#define SIGMA 0
#define BUMP 1
#define BUMP_X 2
#define SWAP 3

/* The total number of plots (used in struct plot) */

#define N_PLOTS 2

/* The gridpoints for plots, usually a fraction of N_GRIDPOINTS */


/*

#ifdef PERIODIC
#define N_GRID_PLOT_PTS_1 (((N_GRIDPOINTS_1) / FACTOR_1))
#undef NO_LAST_POINT
#else
	#ifdef NO_LAST_POINT
	        #define N_GRID_PLOT_PTS_1 (((N_GRIDPOINTS_1-1) / FACTOR_1)) //do not plot the last point
	#else
			#define N_GRID_PLOT_PTS_1 (((N_GRIDPOINTS_1-1) / FACTOR_1) +1)
	#endif
#endif
*/

/* The gridpoints for plots, usually a fraction of N_GRIDPOINTS */

// number of points to be ploted at specific places using adisco_point

#define N_POINTS 3

/* The number PI */
#define PI 3.141592653589793238462643

/* --------------------------------------------------------------*/

/* The size of variables */

#define FLOAT double

/* To print the values of some macros */
/* Warning: only use where OUT_PUT is properly defined */

#define Strg(x) #x
#define PRINT_MACRO_VALUE(x) fprintf(OUT_PUT,"<li> Using %s = %s  </br>\n",#x,Strg(x));
#define GET_MACRO_VALUE(x)   sprintf(macro_value_strg,"%s",Strg(x))


/* -------------------- FUNCTIONS -------------------------------*/

/* different types of plotting routines, could be:
   graph(er-3d), map(le)3, map(le)5, mtv, idl, sv, fsv or gen
   this are changes on the type of file where
   the routine is defined */

#define SDF  /* use to change the plot structure */
#undef SDF
#undef SV
#undef NETCDF
#define ASCHI
#define PYGRAPH

#ifdef SDF
	#include "sdf.h"
	#define ADISCO adisco_sdf_1d
	#define ADISCO_POINT adisco_txt_point_1d
	//#else
	//#define ADISCO adisco_dummy_1d
#endif // SDF

#ifdef ASCHI
	#define ADISCO adisco_aschi_1d
#endif // ASCHI

#ifdef NETCDF
	#define ADISCO adisco_netcdf_1d
#endif // NETCDF

#ifdef SV
	#define ADISCO adisco_fsv_1d   /* #define ADISCO adisco_sv_3d */
#endif // SV

/* -----------------------------------------------------------------*/

/* input output definitions */

#define FILE_INPUT
#define OUT_PUT file_data_ptr

/* used in routines to get data from the web */
#define MAX_ENTRIES 400
#define MAX_ENTRY 200

/* ----------------------- structures -------------------------*/
#define GRID_PAR grid_1d         /* grid parameters */
#define PLOT_PAR plot_1d         /* ploting parameters */
#define INI_PAR w2_ini_par_1d    /* where initial data parameters are stored */
#define FUNCTION_PAR w2_par_1d   /* equation parameters */

/* other structures */
#define INPUT_FUNCTION w2_1d_input_file

/* -----------------------------------------------------------------*/

/* different functions for the integrator, here goes
   most of the physical input */

#define FF w2_eq

/* different arguments for the function FF */

#define FF_ARG struct grid_1d *grid_1d, struct field_array  *fields, struct field_array  *derivs, struct FUNCTION_PAR *function_par
#define FF_ARG_DUMMY struct GRID_PAR *, struct field_array  *, struct field_array  *, struct FUNCTION_PAR *

/* different functions to take derivatives derivQ_1d, derivQ_3_1d, derivD_1d, derivQQ_1d, deriv_strand_third_order_boundaries_sixth_interior, deriv_strand_fourth_order_boundaries_eight_interior, etc */

#ifdef PERIODIC
	// #define DERIV derivD_Per_1d
	// #define DERIV derivQ_Per_1d
	// #define DERIV deriv_6_Per_1d
	#define DERIV deriv_8_Per_1d
#else
	//#define DERIV derivS_1d
	#define DERIV deriv_strand_fourth_order_boundaries_eight_interior_1d
#endif // PERIODIC

#define NON_OPTIMIZED_EIGHT_ORDER_OPERATOR
#undef NON_OPTIMIZED_EIGHT_ORDER_OPERATOR

#define OPTIMIZED_EIGHT_ORDER_OPERATOR
//#undef OPTIMIZED_EIGHT_ORDER_OPERATOR

/* -----------------------------------------------------------------*/

/* different dissipative operators diss_KO_4_D_1d, diss_KO_4_00_D_1d */

#ifdef DISSIPATION
	#define DISS diss_KO4_D_1d
	// #define DISS diss_KO6_Per_1d
	// #define DISS diss_KO8_Per_1d
#endif // DISSIPATION

/* different runge-kutta routines */

//#define RKX rk3
//#define RKX tvd3
#define RKX rk4

/* -----------------------------------------------------------------*/

#endif // FIRST_MACRO_1D_H
