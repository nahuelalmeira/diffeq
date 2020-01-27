/*********************************************************************
*                                                                    *
* This is the Header file of where general structures are defined    * 
*                                                                    *
*                                                                    *
*********************************************************************/

#include "first_macro_1d.h"  /* place where dimensions and float are defined */
#include "globals.h"
#ifdef STRUCTS_1D_H
#else
#define STRUCTS_1D_H

/* ----------------> STRUCTURES <---------------------------------- */


struct field_array {
    /* Value of Field at gridpoint */
    FLOAT *u[N_FIELDS];
    FLOAT time;                 /* time at which the values are taken */
};


struct deriv_array {
    FLOAT *du[N_DERIVS];
    FLOAT time;                 /* time at which the values are taken */
};


struct aux_array {
    FLOAT *u_aux[N_AUX];
    FLOAT time;                 /* time at which the values are taken */
};


struct grid_1d {
  int n_grid_pts; /* Number of gridpoints in coordinate 1*/
  int start_grid;  /* Starting value for grid (ussually = zero) */
  int final_grid;  /* Final value for the grid (ussually = n_gridpts) */

  int n_fields;  /* Number of fields     */
  int data_steps; /* Number of time steps where data is collected */
  int int_steps; /* Number of internal time steps inbetween data collecting */
  int factor_1d_steps;

  FLOAT initial_time; /* Initial time */
  FLOAT final_time; /* Final time */
  FLOAT time_step;
  FLOAT initial_x;
  FLOAT final_x;
  
  int n_rk;        // Runge-Kutta steps
};





#endif
