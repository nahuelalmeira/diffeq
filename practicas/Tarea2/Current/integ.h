/*********************************************************************
*                                                                    *
* integ  -- integrates dy/dt = f(y,t) between two time points        *
*                                                                    *
* Parameters:                                                        *
*   y_b_ptr     -- pointer to the initial data in field_vector struct*
*   t1 (t2)     -- initial (final) integration times                 *
*   data_steps  -- number of steps at which data is saved            *
*   int_steps   -- number of steps between data savings              *
*                  (so total # of steps is data_steps x int_steps    *
*       F       -- pointer to fuction which evaluates the f(y,t)     *
*       rk      -- pointer to runge-kutta (or any other) engine      *
*                                                                    *
* Returns: nothing                                                   *
*                                                                    *
*********************************************************************/

#ifndef __INTEG_H
#define __INTEG_H

#include "first_macro_1d.h"  /* Where global parameters are defined */
#include "structs_1d.h"      /* Where structures are defined */
#include "equation.h"

#ifdef IMEX 

void integ(struct field_array *y_ptr, 
	   struct GRID_PAR *grid_ptr,
	   struct FUNCTION_PAR *equation_par,
	   void (* FF)(FF_ARG_DUMMY),void (* FS)(FF_ARG_DUMMY),void (* FI)(FF_ARG_DUMMY), 
	   void (* RKX)(struct field_array *, 
			struct GRID_PAR *,
			struct FUNCTION_PAR *, 
			void (* )(struct GRID_PAR *,
				  struct field_array *, 
				  struct field_array *,
				  struct FUNCTION_PAR *),
			void (* )(struct GRID_PAR *,
				  struct field_array *, 
				  struct field_array *,
				  struct FUNCTION_PAR *),
			void (* )(struct GRID_PAR *,
				  struct field_array *, 
				  struct field_array *,
				  struct FUNCTION_PAR *)
			)
	   );
	   
#else
	   
void integ(struct field_array *y_ptr,
           struct GRID_PAR *grid_ptr,
           struct FUNCTION_PAR *equation_par,
           void (* FF)(FF_ARG_DUMMY),
           void (* RKX)(struct field_array *,
                        struct GRID_PAR *,
                        struct FUNCTION_PAR *,
                        void (*)(FF_ARG_DUMMY)
                       )
          );

#endif //IMEX




#endif
