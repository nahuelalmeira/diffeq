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

#include "first_macro_1d.h"  /* Where global parameters are defined */
#include "structs_1d.h"      /* Where structures are defined */
#include "derivs_1d.h"       /* Where derivatives functions are defined */
#include "inidat.h"


/* struct field_array *inidat(struct field_array *y_a_ptr) { */

void inidat(struct field_array *y_a_ptr,
	    struct GRID_PAR *grid_1d_ptr,
	    struct INI_PAR *ini_par_ptr){

/* -------> grid parameters <--------------------------*/



    int ni_1 = (*grid_1d_ptr).start_grid;
    int nf_1 = (*grid_1d_ptr).final_grid;

    FLOAT xi = (*grid_1d_ptr).initial_x;
    FLOAT xf = (*grid_1d_ptr).final_x;

/* -------> initial data parameters <------------------*/

#ifdef PERIODIC
 FLOAT twoPIdN1 = 2.*PI*(xf-xi)/(FLOAT)(nf_1-ni_1);
 FLOAT one_dN1 = (xf-xi)/(FLOAT)(nf_1-ni_1);
#else
 FLOAT twoPIdN1 = 2.*PI*(xf-xi)/(FLOAT)(nf_1-ni_1-1);
 FLOAT one_dN1 = (xf-xi)/(FLOAT)(nf_1-ni_1-1);
#endif

  /* Parameters from main */

  FLOAT a0 = (*ini_par_ptr).a0;
  FLOAT k_a_10 = (*ini_par_ptr).k_a_10;
  FLOAT shift_a0 = (*ini_par_ptr).shift_a0;

/* Amplitude sin(k_a1* x + shift_a1) in U1 */

  FLOAT a1 = (*ini_par_ptr).a1;
  FLOAT k_a_11 = (*ini_par_ptr).k_a_11;
  FLOAT shift_a1 = (*ini_par_ptr).shift_a1;

/* Amplitude of cos(k_a0* x + shift_a0) in U0 */

  FLOAT c0 = (*ini_par_ptr).c0;
  FLOAT k_c_10 = (*ini_par_ptr).k_c_10;
  FLOAT shift_c0 = (*ini_par_ptr).shift_c0;

/* Amplitude of cos(k_c1* x + shift_c0) in U1 */
  FLOAT c1 = (*ini_par_ptr).c1;
  FLOAT k_c_11 = (*ini_par_ptr).k_c_11;
  FLOAT shift_c1 = (*ini_par_ptr).shift_c1;

/* Amplitude of b0*exp(-((x-c0_1)^2 + (y-c0_2)^)/sigma_b0) in U0 */
  FLOAT b0 = (*ini_par_ptr).b0;
  FLOAT sigma_b0 = (*ini_par_ptr).sigma_b0;
  FLOAT c0_1 = (*ini_par_ptr).c0_1;


/* Amplitude of exp(cos(k_b1*x)^2/sigma_b1) in U1 */
  FLOAT b1 = (*ini_par_ptr).b1;
  FLOAT sigma_b1 = (*ini_par_ptr).sigma_b1;
  FLOAT c1_1 = (*ini_par_ptr).c1_1;

/* Global wave motion */

  FLOAT v1 = (*ini_par_ptr).v1;

#ifdef DEBUG
  printf("m = %f\n",m);
#endif

/*---------> values for different fields <-------------*/

/* struct field_array fields y; */

/* first the time */

/* y.a.time = initial_time; */

y_a_ptr->time = grid_1d_ptr->initial_time;

{
    register int g_ind1;

    #ifdef SIMPLE_BUMP
        FLOAT x, x_0=0.25, x_1=0.75;
    #endif
    
    #ifdef SQUARE_BUMP
        FLOAT x, x_0=0.25, x_1=0.35, x_2=0.65, x_3=0.75;
    #endif

    FLOAT power = 4.0, power_m_one = power-1.0, powerx2 = 2.0 * power;
    FLOAT b, b_x;

    for (g_ind1 = ni_1; g_ind1 < nf_1; ++g_ind1) {
        x = xi + (FLOAT)g_ind1*one_dN1;


    #ifdef SIMPLE_BUMP // tested ok (using BUMP_TEST)
        if(x >= x_0 && x <= x_1) {
            b   = c0 * pow(x-x_0,power) * pow(x-x_1,power) / pow(0.25,powerx2);
            b_x = c0 * power * (pow(x-x_0,power_m_one) * pow(x-x_1,power) + pow(x-x_0,power) * pow(x-x_1,power_m_one) ) / pow(0.25,powerx2);
      	}
    #endif
/*
    #ifdef SIMPLE_BUMP // tested ok (using BUMP_TEST)
        FLOAT periodo;
        periodo=2*PI*50.0;
        if(x >= x_0 && x <= x_1) {
            b   = sin(periodo*x)*c0 * pow(x-x_0,power) * pow(x-x_1,power) / pow(0.25,powerx2);
            b_x = sin(periodo*x)*c0 * power * (pow(x-x_0,power_m_one) * pow(x-x_1,power) + pow(x-x_0,power) * pow(x-x_1,power_m_one) ) / pow(0.25,powerx2) + periodo*cos(periodo*x)*c0 * pow(x-x_0,power) * pow(x-x_1,power) / pow(0.25,powerx2);
        }
    #endif
*/
#ifdef SQUARE_BUMP // tested ok
    if(x_0 <= x && x <= (x_0+x_1)/2.0){
    b   = c0 * pow(x-x_0,power) * pow(x-x_1,power) / pow((x_1-x_0)/2.,powerx2);
    b_x = c0 * power * (pow(x-x_0,power_m_one) * pow(x-x_1,power) + pow(x-x_0,power) * pow(x-x_1,power_m_one) ) / pow((x_1-x_0)/2.,powerx2);
	}
	else if((x_0+x_1)/2.0 <= x && x <= (x_2+x_3)/2.0){
	b   = 1.0;
	b_x = 0.0;
	}
	else if((x_2+x_3)/2.0<=x && x<=x_3){
	b   = c0 * pow(x-x_2,power) * pow(x-x_3,power) / pow((x_3-x_2)/2.,powerx2);
    b_x = c0 * power * (pow(x-x_2,power_m_one) * pow(x-x_3,power) + pow(x-x_2,power) * pow(x-x_3,power_m_one) ) / pow((x_3-x_2)/2.,powerx2);
	}
#endif

	else{ b = 0.0, b_x = 0.0; }

	globals.auxfields.u_aux[BUMP][g_ind1] = b;
	globals.auxfields.u_aux[BUMP_X][g_ind1] = b_x;

#ifdef BUMP_TEST
	   y_a_ptr->u[U][g_ind1] = b;
	   y_a_ptr->u[V][g_ind1] = b_x;
#endif

//    y_a_ptr->u[U][g_ind1] = a0*sin(k_a_10*x*PI + shift_a1)*exp(-(x-c1_1)*(x-c1_1)*sigma_b1);
//    y_a_ptr->u[V][g_ind1] = a0*k_a_10*PI*cos(k_a_10*x*PI + shift_a1)*exp(-(x-c1_1)*(x-c1_1)*sigma_b1);
	   y_a_ptr->u[U][g_ind1] = b;
	   y_a_ptr->u[V][g_ind1] = 0.0;

    }
  }






  printf("<LI>Inidat finished </br></LI>\n");
}
