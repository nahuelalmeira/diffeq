
/*********************************************************************
*                                                                    *
* This is the Header file of gen.c a general solver for evolutinary  *
* PDE's which uses standard ODE solvers and Spectral derivatives.    *
*                                                                    *
* Basic definitions include the number of grid points, the number of *
* variables.                                                         *
*                                                                    *
*--------------------------------------------------------------------*
*                                                                    *
* Basic structure is of a field:                                     *
*       struc field_array                                            *
* It holds the array pointing to each gridpoint and the time         *
* There will be also an struc bigvector which points to the whole    *
* array structure (a struct field_array)                                          *
*                                                                    *
*--------------------------------------------------------------------*
*                                                                    *
* Functions in different files:                                      *
*                                                                    *
*    intdat -- provides the initial data (it contains physical parm. *
*    odeint -- integrates ordinary diff. equations dy/dt=F(y,t).     *
*    F      -- routine called by odeint which computes the value of  *
*              the right hand side of the above equation. Here are   *
*              physical parameters.                                  *
*    derivs -- computes derivatives using FFT, the relevant subrouti *
*              nes are contained there.                              *
*    adisco -- Send information to the files                         *
*                                                                    *
*                                                                    *
*********************************************************************/

#ifdef GEN_1D_H
#else
#define GEN_1D_H


#include "first_macro_1d.h"  /* Where global parameters are defined */
#include "structs_1d.h"      /* Where structures are defined */
#include "derivs_1d.h"       /* Where derivatives functions are defined */







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

struct w2_par_1d { 
  struct GRID_PAR *grid_ptr;
  FLOAT c;      /* mass of field */
  FLOAT s; 
  FLOAT sigma;   /* dissipation parameter */
  FLOAT R10;     /* Boundary condition at x=0, m = R0*p */
  FLOAT R11;     /* Boundary condition at x=1, p = R1*m */
  FLOAT Omega;   /* Perturbation frequency */
  FLOAT omega;   /* Incomming wave frequency */
  FLOAT lambda;  /* Perturbation strength */


  /* normals to fases, first digit = coord. second = value */

  FLOAT nx_10;
  FLOAT nx_11;

};


/* ---------------------------> FUNCTIONS <------------------------ */


/*********************************************************************
*                                                                    *
* integ  -- integrates dy/dt = f(y,t) between two time points        *
*                                                                    *
* Parameters:                                                        *
*   y_b_ptr     -- pointer to the initial data in field_vector struct*
*       h       -- FLOAT given the time step zice                   *
*   int_steps   -- number of steps between data savings              *
*                  (so total # of steps is data_steps x int_steps    *
*       F       -- pointer to fuction which evaluates the f(y,t)     *
*       RKX     -- pointer to runge-kutta (or any other) engine      *
*                                                                    *
* Returns: nothing                                                   *
*                                                                    *
*********************************************************************/
#ifdef IMAX

extern void integ(struct field_array  *y_b_ptr,  
	   struct GRID_PAR *grid_ptr,
	   struct FUNCTION_PAR *equation_par,
	   void (* FF)(FF_ARG_DUMMY),void (* FS)(FF_ARG_DUMMY),void (* FI)(FF_ARG_DUMMY), 
	   void (* RKX)(struct field_array  *, 
			struct field_array  *,
			struct GRID_PAR *,
			struct FUNCTION_PAR *, 
			void (* )(struct GRID_PAR *,
				  struct field_array  *, 
				  struct field_array  *,
				  struct FUNCTION_PAR *),
			void (* )(struct GRID_PAR *,
				  struct field_array  *, 
				  struct field_array  *,
				  struct FUNCTION_PAR *),
			void (* )(struct GRID_PAR *,
				  struct field_array  *, 
				  struct field_array  *,
				  struct FUNCTION_PAR *)
			)
	   );
#else
extern void integ(struct field_array  *y_b_ptr,
		  struct GRID_PAR *grid_ptr,
		  struct FUNCTION_PAR *function_par,
		  void (* FF)(FF_ARG_DUMMY), 
		  void (* RKX)(struct field_array  *, 
			       struct field_array  *, 
			       struct GRID_PAR *,
			       struct FUNCTION_PAR *, 
			       void (* )(struct GRID_PAR *,
					 struct field_array  *, 
					 struct field_array  *,
					 struct FUNCTION_PAR *)
			       )
		  );
#endif

/*********************************************************************
*                                                                    *
* FF -- evaluates the function f(y,t), lots of physical imput on it  *
*                                                                    *
* Parameters:                                                        *
*       fields -- pointer to field_vector from where to extract (y,t)*
*       derivs -- pointer to field_vector where to put derivatives   *
*       function_par -- pointer to equation parameters               *
* Returns: nothing                                                   *
*                                                                    *
*********************************************************************/

extern void FF(struct GRID_PAR *grid_1d_ptr,
	       struct field_array  *fields, 
	       struct field_array  *derivs, 
	       struct FUNCTION_PAR *function_par
	       );

#ifdef IMAX 
// The complete RHS
extern void FS(struct GRID_PAR *grid_1d_ptr,
	       struct field_array  *fields, 
	       struct field_array  *derivs, 
	       struct FUNCTION_PAR *function_par
	       );
// The implict inversion
extern void FI(struct GRID_PAR *grid_1d_ptr,
	       struct field_array  *fields, 
	       struct field_array  *derivs, 
	       struct FUNCTION_PAR *function_par
	       );
#endif // IMAX








/*********************************************************************
*                                                                    *
* adisco -- sends data to files                                      *
*                                                                    *
* Parameters:                                                        *
*       inst   -- char with instructions ("OPEN", "PUT", or "CLOSE") *
*       fields -- pointer to struct field_array  with data                  *
*                                                                    *
* Returns: execution status                                          *
* Location: file adisco.c (it contains all of them)                  *
*                                                                    *
*********************************************************************/



extern struct PLOT_PAR *adisco_sv_1d(char inst, 
				   struct PLOT_PAR *plot_ptr, 
				   struct GRID_PAR *gri, 
				   struct field_array  *fields);

extern struct PLOT_PAR *adisco_fsv_1d(char inst, 
				   struct PLOT_PAR *plot_ptr, 
				   struct GRID_PAR *gri, 
				   struct field_array  *fields);

extern struct PLOT_PAR *adisco_sdf_1d(char inst, 
				   struct PLOT_PAR *plot_ptr, 
				   struct GRID_PAR *gri, 
				   struct field_array  *fields);



extern struct PLOT_PAR *adisco_txt_point_1d(char inst,
			     struct PLOT_PAR *plot_ptr,
			     struct GRID_PAR *gri,
			     struct field_array  *fields);

extern struct PLOT_PAR *adisco_aschi_1d(char inst,
								 struct PLOT_PAR *plot_ptr,
								 struct GRID_PAR *grid,
								 struct field_array  *fields);


/*********************************************************************
*                                                                    *
* plot_prep   prepares plotting structures to be sent to files       *
*                                                                    *
*********************************************************************/


extern void plot_prep(struct PLOT_PAR *plot_ptr, 
		      struct FUNCTION_PAR *function_par, 
		      struct field_array  *y_ptr
		      );


/* ---------------------------- FF --------------------------------------*/
    
               /* temporal space derivative values */

extern struct field_array_derivs dfields;     
 
extern struct field_array_aux aux_fields;

/********************************************************************/



extern void norm_L2(struct GRID_PAR *grid_1d_ptr,
		struct field_array  *fields_ptr);
		
extern void norm_Energy(struct GRID_PAR *grid_1d_ptr,
		struct field_array  *fields_ptr);

#ifdef LOOPING

extern int p; // for looping
#endif // LOOPING


#endif
