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

#include "first_macro_1d.h"  /* Where global parameters are defined */
#include "structs_1d.h"      /* Where structures are defined */
//#include "derivs_1d.h"       /* Where derivatives functions are defined */
#include "integ.h"
#include "rkc.h"




#ifdef IMAX

void integ(struct field_array *y_b_ptr, 
	   struct GRID_PAR *grid_ptr,
	   struct FUNCTION_PAR *equation_par,
	   void (* FF)(FF_ARG_DUMMY),void (* FS)(FF_ARG_DUMMY),void (* FI)(FF_ARG_DUMMY), 
	   void (* RKX)(struct field_array fields *,
			struct GRID_PAR *,
			struct FUNCTION_PAR *, 
			void (* )(struct GRID_PAR *,
				  struct field_array fields *, 
				  struct field_array fields *,
				  struct FUNCTION_PAR *),
			void (* )(struct GRID_PAR *,
				  struct field_array fields *, 
				  struct field_array fields *,
				  struct FUNCTION_PAR *),
			void (* )(struct GRID_PAR *,
				  struct field_array fields *, 
				  struct field_array fields *,
				  struct FUNCTION_PAR *)
			)
	   ){

printf("revisar integ.c para IMEX!");
exit(0);
/* Total number of components, goten from gen.h */
/*
long int n_total, int_steps;

struct field_array fields *w_ptr;
struct field_array fields w_out;
struct field_array fields *w_out_ptr;

n_total = N_TOTAL;
int_steps = (*grid_ptr).int_steps;
*/
/* get starting values */
/*
w_ptr = y_b_ptr;

w_out_ptr = &w_out;
*/

long int int_steps = (*grid_ptr).int_steps;
 {long int j_inn;
     for (j_inn=1; j_inn<=int_steps; ++j_inn) {/* Take int_steps steps */
	RKX(w_ptr,w_out_ptr,grid_ptr,equation_par,FF,FS,FI); 
	/* swap w with w_out */
	{long int ind;
	for (ind = 0; ind < n_total; ++ind){
	  (*w_ptr).v.u[ind] = (*w_out_ptr).v.u[ind]; 
	}
	}
	(*w_ptr).v.time = (*w_out_ptr).v.time;
/* 	printf("time after RKC in integ = %f",(*w_ptr).v.time); */
      }
 }
}

#else

void integ(struct field_array *y_b_ptr, 
	   struct GRID_PAR *grid_ptr,
	   struct FUNCTION_PAR *equation_par,
	   void (* FF)(FF_ARG_DUMMY), 
	   void (* RKX)(struct field_array *, 
			struct GRID_PAR *,
			struct FUNCTION_PAR *, 
			void (* )(struct GRID_PAR *,
				  struct field_array *, 
				  struct field_array *,
				  struct FUNCTION_PAR *)
			)
	   )
{

	long int int_steps = (*grid_ptr).int_steps;
 	{
		long int j_inn;
     	for (j_inn=1; j_inn<=int_steps; ++j_inn) {/* Take int_steps steps */
			RKX(y_b_ptr, grid_ptr, equation_par, FF); 
			/* 	printf("time after RKC in integ = %f",(*w_ptr).v.time); */
      	}
 	}
}

#endif // IMAX
