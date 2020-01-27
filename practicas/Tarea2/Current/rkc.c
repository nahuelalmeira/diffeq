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

#include "first_macro_1d.h"     /* Where global parameters are defined */
#include "structs_1d.h"         /* Where structures are defined */
//#include "derivs_3d.h"          /* Where derivatives functions are defined */
#include "globals.h"
#include "rkc.h"
#include "assert.h" 


#ifdef IMEX


void imexrkssp3(struct field_array *v_init_ptr, 
	 struct GRID_PAR *grid_ptr,
	 struct FUNCTION_PAR *equation_par,
	 void (* FF)(FF_ARG_DUMMY),void (* FS)(FF_ARG_DUMMY),void (* FI)(FF_ARG_DUMMY) ) 
{ //imexrkssp3

FLOAT x;      /* initial time */   
FLOAT h = (*grid_ptr).time_step; // printf("h=%f\n",h);

int n_fields = (*grid_ptr).n_fields;
long int n_var = (*grid_ptr).n_grid_pts;

FLOAT ae32, ae42, ae43;
FLOAT be2, be3, be4;
FLOAT he1, he2, he3, he4;
FLOAT ai11, ai21, ai22, ai32, ai33, ai41, ai42, ai43, ai44;
FLOAT bi2, bi3, bi4;


//printf("h=%f ",h);

x = v_init_ptr->time;

          ae32 = h;
          ae42 = .25*h;
          ae43 = ae42;
          be2 = h/6.;
          be3 = be2;
          be4 = 2.*h/3.;
          ai11 = 0.24169426078821*h;
          ai21 = -ai11;
          ai22 = ai11;
    	  ai32 = h-ai11;
          ai33 = ai11;
          ai41 = 0.06042356519705*h;
          ai42 = 0.12915286960590*h;
          ai43 = 0.5*h - ai41 - ai42 - ai11;
          ai44 = ai11;
          bi2 = be2;
          bi3 = be3;
          bi4 = be4;
          he1 = 0.0;	
          he2 = 0.0;
          he3 = h;
          he4 = 0.5*h;

// First step, only the invertion

//printf("V0=%f\n",(*v_init_ptr).a.u[V][1]);

// First step

 FI(grid_ptr,v_init_ptr,&globals.v_temp0,equation_par);   // First inversion   
 
 globals.v_temp0.time = x + ai11;

//printf("V1=%1.10f\n",(*v_temp0_ptr).a.u[V][1]);

// U(1) in v_temp0
// No necesitamos F1


// FF(grid_ptr,v_temp0_ptr,dv_tempF1_ptr,equation_par);  //F2

 FS(grid_ptr,&globals.v_temp0,&globals.dv_tempS1,equation_par);  //S(U(1))

//printf("S1=%f\n",(*dv_tempS1_ptr).a.u[V][1]);

#ifdef DEBUG_RKC
printf("after first step in imexrkssp3\n");
#endif

// Second step

	#pragma omp parallel for collapse(2)
	for (int field_ind = 0; field_ind < n_fields; field_ind++){ 
		for (long int ind = 0; ind < n_var; ++ind){
			globals.v_temp0.u[field_ind][ind] = v_init_ptr->u[field_ind][ind] + ai21*globals.dv_tempS1.u[field_ind][ind];
		}
	}


globals.v_temp0.time = x + he2;

 FI(grid_ptr,&globals.v_temp0,&globals.v_temp1,equation_par);   // Second inversion   


//printf("V2=%1.10f\n",(*v_temp1_ptr).a.u[V][1]);

// Third step U(2) in v_temp1 at t=x

 FF(grid_ptr,&globals.v_temp1,&globals.dv_tempF2,equation_par);  // F2 (second step)

 FS(grid_ptr,&globals.v_temp1,&globals.dv_tempS2,equation_par);  // S2...

// U_*(3)


	#pragma omp parallel for collapse(2)
	for (int field_ind = 0; field_ind < n_fields; field_ind++){
		for (long int ind = 0; ind < n_var; ++ind){
			globals.v_temp0.u[field_ind][ind] = (*v_init_ptr).u[field_ind][ind] 
			                                  + ae32*globals.dv_tempF2.u[field_ind][ind]
							  + ai32*globals.dv_tempS2.u[field_ind][ind];
		}
	}
 
#ifdef DEBUG_RKC
printf("after second step in imexrkssp3\n");
#endif

globals.v_temp0.time = x+he3;

 FI(grid_ptr,&globals.v_temp0,&globals.v_temp1,equation_par);   // Third inversion   

//printf("V3=%1.10f\n",(*v_temp1_ptr).a.u[V][1]);

// Fourth step U(3) in v_temp1 at t = x + h

 FF(grid_ptr,&globals.v_temp1,&globals.dv_tempF3,equation_par);  // F3 (third step)

 FS(grid_ptr,&globals.v_temp1,&globals.dv_tempS3,equation_par);  // S3...

// U_*(4)

	#pragma omp parallel for collapse(2)
	for (int field_ind = 0; field_ind < n_fields; field_ind++) { 
		for (long int ind = 0; ind < n_var; ++ind) {
			globals.v_temp0.u[field_ind][ind] = (*v_init_ptr).u[field_ind][ind] 
			                        + ae42*globals.dv_tempF2.u[field_ind][ind]
						+ ae43*globals.dv_tempF3.u[field_ind][ind]
						+ ai41*globals.dv_tempS1.u[field_ind][ind]
						+ ai42*globals.dv_tempS2.u[field_ind][ind]
						+ ai43*globals.dv_tempS3.u[field_ind][ind];
		}
	}

#ifdef DEBUG_RKC
printf("after third step in imexrkssp3\n");
#endif
globals.v_temp0.time = x+he4;

 FI(grid_ptr,&globals.v_temp0,&globals.v_temp1,equation_par);   // Fourth inversion   

//printf("V4=%1.10f\n",(*v_temp1_ptr).a.u[V][1]);

// U(4) in v_temp1 at time t = x + 1/2

 FF(grid_ptr,&globals.v_temp1,&globals.dv_tempF4,equation_par);  // F4 (fourth step)

 FS(grid_ptr,&globals.v_temp1,&globals.dv_tempS1,equation_par);  // S1 in place of S4...

#ifdef DEBUG_RKC
 printf("after calling functions and before adding the final result\n");
#endif

	#pragma omp parallel for collapse(2)
	for (int field_ind = 0; field_ind < n_fields; field_ind++) {
		for (long int ind = 0; ind < n_var; ++ind){
			(*v_init_ptr).u[field_ind][ind] = (*v_init_ptr).u[field_ind][ind]
			                                + be2*(globals.dv_tempF2.u[field_ind][ind] + globals.dv_tempS2.u[field_ind][ind])
							+ be3*(globals.dv_tempF3.u[field_ind][ind] + globals.dv_tempS3.u[field_ind][ind])
							+ be4*(globals.dv_tempF4.u[field_ind][ind] + globals.dv_tempS1.u[field_ind][ind]);
		}
	}


//printf("VF=%1.10f\n",(*v_fina_ptr).a.u[V][1]);

 /* asign returning time */
v_init_ptr->time = x + h;
// printf("time at RKC = %f",x+h); 

#ifdef DEBUG_RKC
printf("exiting imexrkssp3\n");
#endif

//exit(0);
}


// L stable 343 acoording to Asher et al.









/****************** global variables defined above main *****************/

/* union fields dv_tempF2;  */       /* derivative temporal value         */
/* union fields dv_tempF3;  */       /* derivative temporal value         */
/* union fields dv_tempF4;  */       /* derivative temporal value         */
/* union fields dv_tempS1;  */       /* derivative temporal value         */
/* union fields dv_tempS2;  */       /* derivative temporal value         */
/* union fields dv_tempS3;  */       /* derivative temporal value         */

/* union fields dv_sum;   */       /* sum values of deriv               */
/* union fields v_temp0;   */       /* intermediate value of v           */
/* union fields v_temp1;   */       /* intermediate value of v           */

/******************* end global variables *******************************/


// void imexrkL343(union fields *v_init_ptr, 
// 	 union fields *v_fina_ptr, 
// 	 FLOAT h,  
// 	 struct GRID_PAR *grid_ptr,
// 	 struct FUNCTION_PAR *equation_par,
// 	 void (* FF)(FF_ARG_DUMMY),void (* FS)(FF_ARG_DUMMY),void (* FI)(FF_ARG_DUMMY) ) 
// { //imexrkssp3
// 
// FLOAT x;      /* initial time */    
// //FLOAT h, hh;
// 
// long int nvar; /* total number of variables in vector */
// 
// FLOAT ae21, ae31, ae32, ae41, ae42, ae43;
// FLOAT be1, be2, be3, be4;
// FLOAT ai11, ai21, ai22, ai31, ai32, ai33;
// FLOAT bi1, bi2, bi3;
// FLOAT h1, h2, h3;
// 
// 
// printf("Not working YET");
// 
// x = (*v_init_ptr).v.time;
// 
// 
// 
// nvar = N_TOTAL; /* obtained from gen.h */
// 
//  
//           ai11 = 0.4358665215*h;
//           ai21 = (h-ai11)/2.;
// 	  ai22 = ai11;
// 	  ai31 = 3.*ai11*ai11/2./h + 4.*ai11 - 0.25*h;
// 	  ai32 = -ai31 - ai11 + h;
// 	  ai33 = ai11;
//           bi1 = ai31;
//           bi2 = ai32;
//           bi3 = ai11;
// 
//           ae21 = ai11 ;
//           ae31 = 0.3212788860*h;
//           ae32 = 0.3966543747*h;
//           ae41 = -0.105858296*h;
//           ae42 = 0.5529291479*h;
//           ae43 = ae42;
// //           be1 = 0.0;
// // 	  be2 = bi1;
// // 	  be3 = bi2;
// // 	  be4 = bi3;
// 
// 	h1 = ai11*h;
// 	h2 = (1.+ai11)/2.*h;
// 	h3 = h;
// 
// printf("h=%f\n",h);
// 
// // First step, Ke1 = f(u_n) o sea calculo la parte del RHS que no es stiff
// 
// //printf("sin terminar");
// 
// 
//  FF(grid_ptr,v_init_ptr,dv_tempF1_ptr,equation_par);       /* take derivatives at initial time */
// 
//  
// // compute the first intermediate state 
// 
//  {long int ind;
//  for (ind = 0; ind < nvar; ++ind)
//    (*v_temp0_ptr).v.u[ind] = (*v_init_ptr).v.u[ind] + ae21*(*dv_tempF1_ptr).v.u[ind];
//  }
//  
// (*v_temp0_ptr).v.time = x + h1;
// 
//  FI(grid_ptr,v_temp0_ptr,v_temp1_ptr,equation_par); // time is set on the funtion 
// 
//  
// // now we have en v_temp0 the first value and now compute Ki1=g(u(1))
// 
//  FS(grid_ptr,v_temp1_ptr,dv_tempS1_ptr,equation_par);
// 
// 
// // compute the second intermediate step
// 
//  FF(grid_ptr,v_temp1_ptr,dv_tempF2_ptr,equation_par);       /* take derivatives at initial time */
//  
// 
// 
//  {long int ind;
//  for (ind = 0; ind < nvar; ++ind)
//    (*v_temp0_ptr).v.u[ind] = (*v_init_ptr).v.u[ind] + ae31*(*dv_tempF1_ptr).v.u[ind] + ae32*(*dv_tempF2_ptr).v.u[ind] + ai21*(*dv_tempS1_ptr).v.u[ind];
//  }
// 
//  (*v_temp0_ptr).v.time = x + h2;
// 
// 
// 
//  FI(grid_ptr,v_temp0_ptr,v_temp1_ptr,equation_par); // time is set on the funtion 
// 
//  
// // now we have en v_temp0 the first value and now compute Ki1=g(u(1))
// 
//  FS(grid_ptr,v_temp1_ptr,dv_tempS2_ptr,equation_par);
// 
// 
// // compute the third intermediate step
// 
//  FF(grid_ptr,v_temp1_ptr,dv_tempF3_ptr,equation_par);       /* take derivatives at initial time */
//  
//  {long int ind;
//  for (ind = 0; ind < nvar; ++ind)
//    (*v_temp0_ptr).v.u[ind] = (*v_init_ptr).v.u[ind] + ae41*(*dv_tempF1_ptr).v.u[ind] + ae42*(*dv_tempF2_ptr).v.u[ind] + ae43*(*dv_tempF3_ptr).v.u[ind] + ai31*(*dv_tempS1_ptr).v.u[ind] + ai32*(*dv_tempS2_ptr).v.u[ind];
//  }
// 
//  (*v_temp0_ptr).v.time = x + h3;
// 
// // compute the last not stiff RHS
// 
//  FF(grid_ptr,v_temp1_ptr,dv_tempF4_ptr,equation_par);       /* take derivatives at initial time */
// 
//  {long int ind;
//  for (ind = 0; ind < nvar; ++ind){
//    (*v_fina_ptr).v.u[ind] = (*v_init_ptr).v.u[ind] + bi1*((*dv_tempS1_ptr).v.u[ind] + (*dv_tempF2_ptr).v.u[ind]) + bi2*((*dv_tempS2_ptr).v.u[ind] + (*dv_tempF3_ptr).v.u[ind]) + bi3*((*dv_tempS3_ptr).v.u[ind] + (*dv_tempF4_ptr).v.u[ind]);
//  }}
// 
// 
// 
//  /* asign returning time */
// (*v_fina_ptr).v.time = x + h;
// /*  printf("time at RKC = %f",x+h); */
// 
// 
// 
// 
// }




#else

void rk4(   struct field_array *v_init_ptr, struct GRID_PAR *grid_ptr, struct FUNCTION_PAR *equation_par, void (*FF) (FF_ARG_DUMMY))
{
    FLOAT h = (*grid_ptr).time_step; // printf("h=%f\n",h);
    FLOAT x;                    /* initial time */
    FLOAT xh;                   /* initial time + hh */
    FLOAT hh;                   /* hh half the h */
    FLOAT h6;                   /* one sixth of h */
    int n_fields = grid_ptr->n_fields;
    long int n_var = (*grid_ptr).n_grid_pts;

/****************** global variables defined above main *****************/


    /* struct field_array dv_temp;  *//* derivative temporal value         */
    /* struct field_array dv_sum;   *//* sum values of deriv               */
    /* struct field_array v_temp;   *//* intermediate value of v           */

    /* struct field_array *dv_temp_ptr = &dv_temp;  *//* derivative temporal value        */
    /* struct field_array *dv_sum_ptr = &dv_sum;    *//* sum          values of deriv     */
    /* struct field_array *v_temp_ptr = &v_temp;    *//* intermediate value of v          */


/******************* end global variables *******************************/

    x = v_init_ptr->time;
    hh = h * (0.5);
    h6 = h / (6.0);
    xh = x + hh;


    /* take derivatives at initial time */
    grid_ptr->n_rk = 0;
    FF(grid_ptr, v_init_ptr, &globals.dv_temp, equation_par);

#ifdef DEBUG
    printf("after first function evaluation in RKC \n");
#endif                          /* DEBUG */

    /* put derivative in dev_sum = (k1 + 2k2 + 2k3 + k4) */
    /* dev_sum = k1 */
    for (int field_ind = 0; field_ind < n_fields; field_ind++) {
        for (long ind = 0; ind < n_var; ++ind) {
            globals.dv_sum.u[field_ind][ind] = globals.dv_temp.u[field_ind][ind];
        }
    }

    // return(0); 
    // printf("me equivoque");

    /* asign values to v_temp (v_temp = v_init + hh*dv_init) */
    for (int field_ind = 0; field_ind < n_fields; field_ind++) {
        for (long ind = 0; ind < n_var; ++ind) {
            globals.v_temp.u[field_ind][ind] = v_init_ptr->u[field_ind][ind] + globals.dv_temp.u[field_ind][ind] * hh;
                      //assert(globals.v_temp.u[field_ind][ind] < 100000.0);
        }
    }


    /* asign time to v_temp */
    globals.v_temp.time = xh;

    /* take again derivative at time xh = t + h/2 and put them in dv_temp */
    grid_ptr->n_rk = 1;
    FF(grid_ptr, &globals.v_temp, &globals.dv_temp, equation_par);

    /* put derivative in dev_sum = (k1 + 2k2 + 2k3 + k4) */
    /* dev_sum = k1 + 2k2 */
    for (int field_ind = 0; field_ind < n_fields; field_ind++) {
        for (long ind = 0; ind < n_var; ++ind) {
            globals.dv_sum.u[field_ind][ind] += 2 * globals.dv_temp.u[field_ind][ind];
        }
    }

    /* assign values to v_temp = v_init + (h/2)*dv_temp  */
    for (int field_ind = 0; field_ind < n_fields; field_ind++) {
        for (long ind = 0; ind < n_var; ++ind) {
            globals.v_temp.u[field_ind][ind] = v_init_ptr->u[field_ind][ind] + globals.dv_temp.u[field_ind][ind] * hh;
                      // assert(globals.v_temp.u[field_ind][ind] < 100000.0);
        }
    }

    /* time does not need to be updated */

    /* takes derivatives at x+h/2, v_temp (new) and put them in dv_med */
    grid_ptr->n_rk = 2;
    FF(grid_ptr, &globals.v_temp, &globals.dv_temp, equation_par);

    /* put derivative in dev_sum = (k1 + 2k2 + 2k3 + k4) */
    /* dev_sum = k1 + 2k2 + 2k3 */
    for (int field_ind = 0; field_ind < n_fields; field_ind++) {
        for (long ind = 0; ind < n_var; ++ind) {
            globals.dv_sum.u[field_ind][ind] += 2 * globals.dv_temp.u[field_ind][ind];
        }
    }

    /* asign values to v_temp = v_init + h*dv_temp */
    for (int field_ind = 0; field_ind < n_fields; field_ind++) {
        for (long ind = 0; ind < n_var; ++ind) {
            globals.v_temp.u[field_ind][ind] = v_init_ptr->u[field_ind][ind] + globals.dv_temp.u[field_ind][ind] * h;
                      //  assert(globals.v_temp.u[field_ind][ind] < 100000.0);
        }
    }

    /* asigns new time to v_temp = x + h */
    globals.v_temp.time = x + h;

    /* takes derivatives at x+h with new v_temp and puts them in dv_temp */
    grid_ptr->n_rk = 3;
    FF(grid_ptr, &globals.v_temp, &globals.dv_temp, equation_par);

    /* asign value to v_fina = v_init + h6*(dv_temp + dv_temp) */
    for (int field_ind = 0; field_ind < n_fields; field_ind++) {
        for (long ind = 0; ind < n_var; ++ind) {
            v_init_ptr->u[field_ind][ind] += h6 * (globals.dv_sum.u[field_ind][ind] + globals.dv_temp.u[field_ind][ind]);
                       // assert(v_init_ptr->u[field_ind][ind] < 100000.0);
        }
    }

    /* asign returning time */
    v_init_ptr->time = x + h;

    // printf("time at RKC = %f",x+h);
}


void rk3(   struct field_array *v_init_ptr, struct GRID_PAR *grid_ptr, struct FUNCTION_PAR *equation_par, void (*FF) (FF_ARG_DUMMY))
{
	
/***********************************************************************************
 *                                                                                 *
 *   Third order Runge Kutta                                                       *
 *                                                                                 *
 **********************************************************************************/
    FLOAT h = (*grid_ptr).time_step; // printf("h=%f\n",h);
    FLOAT x;                    /* initial time */
    FLOAT xh;                   /* initial time + hh */
    FLOAT hh;                   /* hh half the h */
    FLOAT h9;                   /* one ninth of h */
    FLOAT h34;

    int n_fields = grid_ptr->n_fields;
    long int n_var = (*grid_ptr).n_grid_pts;

/****************** global variables defined above main *****************/


    /* struct field_array dv_temp;  *//* derivative temporal value         */
    /* struct field_array dv_sum;   *//* sum values of deriv               */
    /* struct field_array v_temp;   *//* intermediate value of v           */

    /* struct field_array *dv_temp_ptr = &dv_temp;  *//* derivative temporal value        */
    /* struct field_array *dv_sum_ptr = &dv_sum;    *//* sum          values of deriv     */
    /* struct field_array *v_temp_ptr = &v_temp;    *//* intermediate value of v          */


/******************* end global variables *******************************/

    x = v_init_ptr->time;
    hh = h * (0.5);
    h9 = h / (9.0);
    h34 = h * 3.0 / 4.0;
    xh = x + hh;

    /* take derivatives at initial time */
    grid_ptr->n_rk = 0;
    FF(grid_ptr, v_init_ptr, &globals.dv_temp, equation_par);

#ifdef DEBUG
    printf("after first function evaluation in RKC \n");
#endif

    /* put derivative in dev_sum = (2*k1 + 3*k2 + 4*k3) */
    /* dev_sum = 2*k1 */
    for (int field_ind = 0; field_ind < n_fields; field_ind++) {
        for (long ind = 0; ind < n_var; ++ind) {
            globals.dv_sum.u[field_ind][ind] = 2.0 * globals.dv_temp.u[field_ind][ind];
        }
    }

    /* asign values to v_temp (v_temp = v_init + hh*dv_init) */
    for (int field_ind = 0; field_ind < n_fields; field_ind++) {
        for (long ind = 0; ind < n_var; ++ind) {
            globals.v_temp.u[field_ind][ind] = v_init_ptr->u[field_ind][ind] + hh * globals.dv_temp.u[field_ind][ind];
        }
    }

    /* asign time to v_temp */
    globals.v_temp.time = xh;

    /* take again derivative at time xh = t + h/2 and put them in dv_temp */
    grid_ptr->n_rk = 1;
    FF(grid_ptr, &globals.v_temp, &globals.dv_temp, equation_par);

    /* put derivative in dev_sum = (2*k1 + 3*k2 + 4*k3) */
    /* dev_sum = 2*k1 + 3*k2 */
    for (int field_ind = 0; field_ind < n_fields; field_ind++) {
        for (long ind = 0; ind < n_var; ++ind) {
            globals.dv_sum.u[field_ind][ind] += 3.0 * globals.dv_temp.u[field_ind][ind];
        }
    }

    /* asign values to v_temp = v_init + (3*h/4)*dv_temp  */
    for (int field_ind = 0; field_ind < n_fields; field_ind++) {
        for (long ind = 0; ind < n_var; ++ind) {
            globals.v_temp.u[field_ind][ind] = v_init_ptr->u[field_ind][ind] + h34 * globals.dv_temp.u[field_ind][ind];
        }
    }

    /* assign time to v_temp */
    globals.v_temp.time = x + h34;

    /* takes derivatives at x+3*h/4, v_temp (new) and put them in dv_med */
    grid_ptr->n_rk = 2;
    FF(grid_ptr, &globals.v_temp, &globals.dv_temp, equation_par);

    /* assign value to v_fina = v_init + h9*(dv_sum + 4*dv_temp) */
    for (int field_ind = 0; field_ind < n_fields; field_ind++) {
        for (long ind = 0; ind < n_var; ++ind) {
            v_init_ptr->u[field_ind][ind] += h9 * (globals.dv_sum.u[field_ind][ind] + 4.0 * globals.dv_temp.u[field_ind][ind]);
        }
    }

    /* asign returning time */
    v_init_ptr->time = x + h;

    // printf("time at RKC = %f",x+h);
}


void rk3lo9(struct field_array *v_init_ptr, struct GRID_PAR *grid_ptr, struct FUNCTION_PAR *equation_par, void (*FF) (FF_ARG_DUMMY))
{


/***********************************************************************************
 *                                                                                 *
 *   Low storage third order Runge Kutta, according to                             *
 *   SIAM J. SCI. COMPUT. Vol. 30, No. 4, pp. 2113–2136                            *
 *   HIGHLY EFFICIENT STRONG STABILITY-PRESERVING RUNGE–KUTTA METHODS              *
 *   WITH LOW-STORAGE                                                              *
 *   IMPLEMENTATIONS                   DAVID I. KETCHESON                          *
 *                                                                                 *
 **********************************************************************************/

    FLOAT h = (*grid_ptr).time_step; // printf("h=%f\n",h);
    FLOAT x;                    /* initial time */
    int n = 3, s = n * n, r = s - n, i;
    FLOAT dt = h / (FLOAT) (r);

    int n_fields = grid_ptr->n_fields;
    long int n_var = (*grid_ptr).n_grid_pts;

    grid_ptr->n_rk = 0;

    /****************** global variables defined above main *****************/


    /* struct field_array dv_temp;  *//* derivative temporal value         */
    /* struct field_array dv_sum;   *//* sum values of deriv               */
    /* struct field_array v_temp;   *//* intermediate value of v           */

    /* struct field_array *dv_temp_ptr = &dv_temp;  *//* derivative temporal value        */
    /* struct field_array *dv_sum_ptr = &dv_sum;    *//* sum          values of deriv     */
    /* struct field_array *v_temp_ptr = &v_temp;    *//* intermediate value of v          */


    /******************* end global variables *******************************/


    // x = (*v_init_ptr).time;

    for (i = 1; i <= (n - 1) * (n - 2) / 2; i++) {

        //q1 = q1 + dt*F(q1)/r; 
        FF(grid_ptr, v_init_ptr, &globals.dv_temp, equation_par);    /* take derivatives at initial time */
        grid_ptr->n_rk++;

#ifdef DEBUG
        printf("after first function evaluation in rk3lo9 \n");
#endif

        for (int field_ind = 0; field_ind < n_fields; field_ind++) {
            for (long ind = 0; ind < n_var; ++ind) {
                v_init_ptr->u[field_ind][ind] = v_init_ptr->u[field_ind][ind] + dt * globals.dv_temp.u[field_ind][ind];
            }
        }
        v_init_ptr->time = v_init_ptr->time + dt;
    }

    // q2=q1;


    /* asign values to v_temp (v_temp = v_init) */

    for (int field_ind = 0; field_ind < n_fields; field_ind++) {
        for (long ind = 0; ind < n_var; ++ind) {
            globals.v_temp.u[field_ind][ind] = v_init_ptr->u[field_ind][ind];
            globals.v_temp.time = v_init_ptr->time;
        }
    }


    // for i=(n-1)*(n-2)/2+1:n*(n+1)/2-1
    // q1 = q1 + dt*F(q1)/r; 
    // end
    for (i = (n - 1) * (n - 2) / 2 + 1; i <= n * (n + 1) / 2 - 1; i++) {

        FF(grid_ptr, v_init_ptr, &globals.dv_temp, equation_par);
        grid_ptr->n_rk++;

        for (int field_ind = 0; field_ind < n_fields; field_ind++) {
            for (long ind = 0; ind < n_var; ++ind) {
                v_init_ptr->u[field_ind][ind] = v_init_ptr->u[field_ind][ind] + dt * globals.dv_temp.u[field_ind][ind];
            }
        }
        v_init_ptr->time = v_init_ptr->time + dt;
    }

    // q1 = ( n*q2 + (n-1)*(q1 + dt*F(q1)/r) ) / (2*n-1);
    FF(grid_ptr, v_init_ptr, &globals.dv_temp, equation_par);
    grid_ptr->n_rk++;
    for (int field_ind = 0; field_ind < n_fields; field_ind++) {
        for (long ind = 0; ind < n_var; ++ind) {
            v_init_ptr->u[field_ind][ind] = ( (FLOAT) (n) * globals.v_temp.u[field_ind][ind]
                                            + (FLOAT) (n - 1) * ( v_init_ptr->u[field_ind][ind]
                                                                + dt * globals.dv_temp.u[field_ind][ind]
                                                                )
                                            ) / (FLOAT) (2 * n - 1);
        }
    }
    v_init_ptr->time = ((FLOAT) (n) * globals.v_temp.time + (FLOAT) (n - 1) * (v_init_ptr->time + dt)) / (FLOAT) (2 * n - 1);

    // for i=n*(n+1)/2+1:s
    // q1 = q1 + dt*F(q1)/r; 
    // end
    for (i = n * (n + 1) / 2 + 1; i <= s; i++) {

        FF(grid_ptr, v_init_ptr, &globals.dv_temp, equation_par);
        grid_ptr->n_rk++;
        for (int field_ind = 0; field_ind < n_fields; field_ind++) {
            for (long ind = 0; ind < n_var; ++ind) {
                v_init_ptr->u[field_ind][ind] = v_init_ptr->u[field_ind][ind] + dt * globals.dv_temp.u[field_ind][ind];
            }
        }
        v_init_ptr->time = v_init_ptr->time + dt;
    }

    //printf("time at RKC = %f, x+h=%f",v_init_ptr->time, x+h);  
    /* asign returning time */
    //v_init_ptr->time = x + h;
    /*  printf("time at RKC = %f",x+h); */
}



void rk2(   struct field_array *v_init_ptr, struct GRID_PAR *grid_ptr, struct FUNCTION_PAR *equation_par, void (*FF) (FF_ARG_DUMMY))
{
	
/***********************************************************************************
 *                                                                                 *
 *   Second order Runge Kutta                                                      *
 *   Mid point method                                                              *
 *                                                                                 *
 **********************************************************************************/
    FLOAT h = (*grid_ptr).time_step; // printf("h=%f\n",h);
    FLOAT hh;
    FLOAT x;                    /* initial time */
    FLOAT xh;                   /* initial time + hh */
    FLOAT a = 0.5;              // parameter to fix method

    int n_fields = grid_ptr->n_fields;
    long int n_var = (*grid_ptr).n_grid_pts;

/****************** global variables defined above main *****************/


    /* struct field_array dv_temp;  *//* derivative temporal value         */
    /* struct field_array v_temp;   *//* intermediate value of v           */

    /* struct field_array *dv_temp_ptr = &dv_temp;  *//* derivative temporal value        */
    /* struct field_array *v_temp_ptr = &v_temp;    *//* intermediate value of v          */


/******************* end global variables *******************************/

    x = v_init_ptr->time;
    hh = h * a;
    xh = x + h * a;

    /* take derivatives at initial time */
    grid_ptr->n_rk = 0;
    FF(grid_ptr, v_init_ptr, &globals.dv_temp, equation_par);

#ifdef DEBUG
    printf("after first function evaluation in RKC \n");
#endif


    /* asign values to v_temp (v_temp = v_init + hh*dv_init) */
    for (int field_ind = 0; field_ind < n_fields; field_ind++) {
        for (long ind = 0; ind < n_var; ++ind) {
            globals.v_temp.u[field_ind][ind] = v_init_ptr->u[field_ind][ind] + hh * globals.dv_temp.u[field_ind][ind];
        }
    }

    /* asign time to v_temp */
    globals.v_temp.time = xh;

    /* take again derivative at time xh = t + h/2 and put them in dv_temp */
    grid_ptr->n_rk = 1;
    FF(grid_ptr, &globals.v_temp, &globals.dv_temp, equation_par);


    /* assign value to v_fina = v_init + h*(dv_temp) */
    for (int field_ind = 0; field_ind < n_fields; field_ind++) {
        for (long ind = 0; ind < n_var; ++ind) {
            v_init_ptr->u[field_ind][ind] += h * globals.dv_temp.u[field_ind][ind];
        }
    }

    /* asign returning time */
    v_init_ptr->time = x + h;

    // printf("time at RKC = %f",x+h);
}




void euler(   struct field_array *v_init_ptr, struct GRID_PAR *grid_ptr, struct FUNCTION_PAR *equation_par, void (*FF) (FF_ARG_DUMMY))
{
	
/***********************************************************************************
 *                                                                                 *
 *   Euler method                                                                  *
 *                                                                                 *
 **********************************************************************************/
 
    FLOAT h = (*grid_ptr).time_step; // printf("h=%f\n",h);
    FLOAT x;                    /* initial time */
    FLOAT xh;                   /* initial time + hh */

    int n_fields = grid_ptr->n_fields;
    long int n_var = (*grid_ptr).n_grid_pts;

/****************** global variables defined above main *****************/


    /* struct field_array dv_temp;  *//* derivative temporal value         */

    /* struct field_array *dv_temp_ptr = &dv_temp;  *//* derivative temporal value        */


/******************* end global variables *******************************/

    x = v_init_ptr->time;
    xh = x + h;

    /* take derivatives at initial time */
    grid_ptr->n_rk = 0;
    FF(grid_ptr, v_init_ptr, &globals.dv_temp, equation_par);

#ifdef DEBUG
    printf("after first function evaluation in RKC \n");
#endif

    /* assign value to v_fina = v_init + h9*(dv_sum + 4*dv_temp) */
    for (int field_ind = 0; field_ind < n_fields; field_ind++) {
        for (long ind = 0; ind < n_var; ++ind) {
            v_init_ptr->u[field_ind][ind] += h * globals.dv_temp.u[field_ind][ind];
        }
    }

    /* asign returning time */
    v_init_ptr->time = x + h;

    // printf("time at RKC = %f",x+h);
}


#endif //IMEX
