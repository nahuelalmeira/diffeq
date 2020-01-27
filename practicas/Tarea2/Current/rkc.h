/*********************************************************************
*                                                                    *
*  rkc -- integrates    dy/dt = f(y,t)     using runge-kutta         *
*                                                                    *
* Parameters:                                                        *
*      v_init_ptr     -- pointer to where the y initial is           *
*                        and where the y final goes                  *
*      h              -- time step                                   *
*      f              -- pointer to function which computes f(y,t)   *
*                                                                    *
* Returns: nothing                                                   *
*                                                                    *
* WARNING uses global variables for intermediate values              *
*********************************************************************/
#ifndef __RKC_H
#define __RKC_H

#include "first_macro_1d.h"     /* Where global parameters are defined */
#include "structs_1d.h"         /* Where structures are defined */
#include "rkc.h"
#include "equation.h"

#ifdef IMEX 

void imexrkssp3(struct field_array *v_init_ptr, 
				struct GRID_PAR *grid_ptr, 
				struct FUNCTION_PAR *equation_par, 
				void (* FF)(FF_ARG_DUMMY), 
				void (* FS)(FF_ARG_DUMMY),
				void (* FI)(FF_ARG_DUMMY)
				);
void imexrkL343(struct field_array  *v_init_ptr, 
	 struct GRID_PAR *grid_ptr,
	 struct FUNCTION_PAR *equation_par,
	 void (* FF)(FF_ARG_DUMMY),
	 void (* FS)(FF_ARG_DUMMY),
	 void (* FI)(FF_ARG_DUMMY) );    

#else

void rk4(   struct field_array *v_init_ptr, struct GRID_PAR *grid_ptr, struct FUNCTION_PAR *equation_par, void (*FF) (FF_ARG_DUMMY));
void rk3(   struct field_array *v_init_ptr, struct GRID_PAR *grid_ptr, struct FUNCTION_PAR *equation_par, void (*FF) (FF_ARG_DUMMY));
void rk2(   struct field_array *v_init_ptr, struct GRID_PAR *grid_ptr, struct FUNCTION_PAR *equation_par, void (*FF) (FF_ARG_DUMMY));
void euler(   struct field_array *v_init_ptr, struct GRID_PAR *grid_ptr, struct FUNCTION_PAR *equation_par, void (*FF) (FF_ARG_DUMMY));
void rk3lo9(struct field_array *v_init_ptr, struct GRID_PAR *grid_ptr, struct FUNCTION_PAR *equation_par, void (*FF) (FF_ARG_DUMMY));

/*********************************************************************
*                                                                    *
*  tvd3 -- integrates    dy/dt = f(y,t)     total variation          *
*                                    diminishing runge-kutta         *
*                                          third order               *
* Parameters:                                                        *
*      v_init_ptr     -- pointer to where the y initial is           *
*      v_fina_ptr     -- pointer to where the y final is             *
*      h              -- time step                                   *
*      function_par   -- pointer to parameters to pass to FF         *
*      FF              -- pointer to function which computes f(y,t)  *
*                                                                    *
* Returns: nothing                                                   *
*                                                                    *
*********************************************************************/


/************************************************************************/
/* The scheme is:                                                       */
/* u^1 = u^n + h*F(u^n,t)                                               */
/* u^2 = 3u^n/4 + u^1/4 + h/4*F(u^1,t+h)                                */
/* u^(n+1) = u^n/3 + 2*u^2/3 * 2h/3*F(u^2,t+h/2)                        */
/*                                                                      */
/* ver paper de Luis del 2008                                           */
/************************************************************************/

void tvd3(struct field_array  *v_init_ptr, struct GRID_PAR *grid_1d_ptr, struct FUNCTION_PAR *function_par,void (* FF)(FF_ARG_DUMMY));



#endif // IMEX
#endif
