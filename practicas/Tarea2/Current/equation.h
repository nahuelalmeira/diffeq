#ifndef __CPU_RICCI_H
#define __CPU_RICCI_H

#include "first_macro_1d.h"     /* Where global parameters are defined */
#include "structs_1d.h"         /* Where structures are defined */




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



#ifdef IMEX

void FF(struct GRID_PAR *grid_ptr, struct field_array *fields_ptr, struct field_array *derivs_ptr, struct FUNCTION_PAR *function_par);
void FS(struct GRID_PAR *grid_ptr, struct field_array *fields_ptr, struct field_array *derivs_ptr, struct FUNCTION_PAR *function_par);
void FI(struct GRID_PAR *grid_ptr, struct field_array *fields_ptr, struct field_array *derivs_ptr, struct FUNCTION_PAR *function_par);

#else

void FF(struct GRID_PAR *grid_ptr, struct field_array *fields_ptr, struct field_array *derivs_ptr, struct FUNCTION_PAR *function_par);

#endif // IMEX







#endif /* __CPU_RICCI_H */
