/*********************************************************************
 *                                                                    *
 * inidat -- provides initial data                                    *
 *                                                                    *
 * Parameters:                                                        *
 * y_a_ptr        -- pointer where to write the initial data          *
 * initial_time   -- initial time                                     *
 *                                                                    *
 * Returns: pointer to field_array where data was writen              *
 *                                                                    *
 *********************************************************************/

#ifndef __INIDAT_H
#define __INIDAT_H

#include "first_macro_1d.h"
#include "structs_1d.h"         /* Where structures are defined */

struct w2_ini_par_1d {
  struct FUNCTION_PAR *function_par_ptr;
  struct GRID_PAR *grid_ptr;
/* Initial data is given by:                                    */
/*       U[0] = PHI = (a0*sin(k_a_10 x + k_a_20 y + shift_a0)   */
/*       + c_0*cos(k_c_10 x + k_c_20 y + shift_c0))*            */
/*       b_0*exp(-sigma_b0*((x-c0_1)^2+(y-c0_2)^2)             */
/* U[1] = dPHI/dt = (a1*sin(k_a_11 x + k_a_21 y + shift_a1)     */
/*       + c_1*cos(k_c_11 x + k_c_21 y + shift_c1))*            */
/*       b_1*exp(-sigma_b1*((x-c1_1)^2+(y-c1_2)^2+)             */
/*       + v1*U[0]_x + v2*U[0]_y                               */

  FLOAT a0;                
  FLOAT k_a_10; 
  FLOAT shift_a0;         

  FLOAT a1;              
  FLOAT k_a_11;   
  FLOAT shift_a1;     

  FLOAT c0;             
  FLOAT k_c_10;
  FLOAT shift_c0;     

  FLOAT c1;            
  FLOAT k_c_11;  
  FLOAT shift_c1;       
  

  FLOAT b0;           
  FLOAT sigma_b0;
  FLOAT c0_1;

  FLOAT b1;               
  FLOAT sigma_b1;
  FLOAT c1_1;

    FLOAT v1;
};

	
void inidat(struct field_array *y_a_ptr, struct GRID_PAR *grid_ptr, struct INI_PAR *ini_par_ptr);

#endif
