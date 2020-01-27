
/*********************************************************************
*                                                                    *
* F -- evaluates the function f(y,t), lots of physical imput on it   *
*                                                                    *
* Parameters:                                                        *
*       fields -- pointer to field_vector from where to extract (y,t)*
*       derivs -- pointer to field_vector where to put derivatives   *
*     wave_par -- pointer to parameter struct                        *
*                                                                    *
* Returns: nothing                                                   *
*                                                                    *
*********************************************************************/
#include "first_macro_1d.h"  /* Where global parameters are defined */
#include "structs_1d.h"      /* Where structures are defined */
#include "derivs_1d.h"       /* Where derivatives functions are defined */
#include "equation.h"
//#include "gen_1d.h"


/***********************************************************************/

void w2_eq(struct GRID_PAR *grid_1d_ptr,
		       struct field_array *fields_ptr,
		       struct field_array *derivs_ptr,
		       struct FUNCTION_PAR *function_par) 
{

    int ni_1 = (*grid_1d_ptr).start_grid;
    int nf_1 = (*grid_1d_ptr).final_grid;

    FLOAT xi = (*grid_1d_ptr).initial_x;
    FLOAT xf = (*grid_1d_ptr).final_x;
    FLOAT dt = (*grid_1d_ptr).time_step;

    FLOAT cc = (*function_par).c;
    FLOAT s = (*function_par).s;

    FLOAT sigma = (*function_par).sigma;
    FLOAT Omega = (*function_par).Omega;
    Omega = 2.0 * PI * Omega;
    FLOAT omega = (*function_par).omega;
    omega = 2.0 * PI * omega;
    FLOAT lambda = (*function_par).lambda;
    //	printf("Omega = %f, omega = %f, lambda = %f\n ", Omega, omega, lambda);
    FLOAT u, v, a, a2, a_t, a_x;
    /* normals */

    FLOAT nx_10 = (*function_par).nx_10;
    FLOAT nx_11 = (*function_par).nx_11;

    FLOAT nx;

    FLOAT factor, L, LL, dU, dV;

    #ifndef PERIODIC

        char macro_value_strg[100];

        GET_MACRO_VALUE(DERIV);
        if (strcmp(macro_value_strg,"derivS_1d")==0) {
          factor = 2.0;
        }
        else if (strcmp(macro_value_strg,"derivQ_1d")==0) {
          factor = 48.0/17.0;
        }
        else if (strcmp(macro_value_strg,"derivQ_3_1d")==0) {
          factor = 11.0/3.0;
        }
        else if (strcmp(macro_value_strg,"deriv_strand_third_order_boundaries_sixth_interior_1d")==0) {
          factor = 43200.0/13649.0;
        }
        else if (strcmp(macro_value_strg,"deriv_strand_fourth_order_boundaries_eight_interior_1d")==0) {
          factor = 5080320.0/1498139.0;
        }
        else {
          factor=2.0; printf("check factor por penalty!!!%s!!!",macro_value_strg);
        }

        L=(FLOAT)(nf_1-ni_1-1)/(xf-xi);
        LL=factor*L;


    #endif

    /* first the time */

    derivs_ptr->time = fields_ptr->time;

    /* inner points */  /* field 0 u1_t = u2 */

    /* take derivatives */
    #ifdef STRICT
        {
            register int grid_ind1;
            #pragma omp parallel for
            for (grid_ind1 = ni_1; grid_ind1 < nf_1; ++grid_ind1) {
                a = 1.0 + lambda * globals.auxfields.u_aux[BUMP][grid_ind1] * cos(PI * Omega * fields_ptr->time);
                a2 = a * a;

                globals.auxfields.u_aux[SWAP][grid_ind1] = fields->u[U][grid_ind1] * a2;
            }
        }

        DERIV(grid_1d_ptr, globals.auxfields.u_aux[SWAP][0],
              globals.dfields.du[U_X][0]);

        {
            register int grid_ind1;
            #pragma omp parallel for
            for (grid_ind1 = ni_1; grid_ind1 < nf_1; ++grid_ind1) {
	              a = 1.0 + lambda * globals.auxfields.u_aux[BUMP][grid_ind1] * cos(Omega * fields_ptr->time);
	              a2 = a * a;

	              globals.auxfields.u_aux[SWAP][grid_ind1] = fields->u[V][grid_ind1] * a2;
            }
        }

        DERIV(grid_1d_ptr, globals.auxfields.u_aux[SWAP], globals.dfields.du[V_X]);

    #else
        DERIV(grid_1d_ptr, fields_ptr->u[U], globals.dfields.du[U_X]);
        DERIV(grid_1d_ptr, fields_ptr->u[V], globals.dfields.du[V_X]);
    #endif // STRICT

    #ifdef DISSIPATION
        DISS(grid_1d_ptr, fields_ptr->u[U], globals.dfields.du[DISS_U]);
        DISS(grid_1d_ptr, fields_ptr->u[V], globals.dfields.du[DISS_V]);
    #endif // DISSIPATION

    /* inner points */
    {
        register int grid_ind1;
        #pragma omp parallel for
        for (grid_ind1 = ni_1; grid_ind1 < nf_1; ++grid_ind1){
            u = fields_ptr->u[U][grid_ind1];
            v = fields_ptr->u[V][grid_ind1];
            a = 1.0 + lambda * globals.auxfields.u_aux[BUMP][grid_ind1] * cos(Omega * fields_ptr->time);
            a_t =   - lambda * Omega * globals.auxfields.u_aux[BUMP][grid_ind1] * sin(Omega * fields_ptr->time);
            a_x =     lambda * globals.auxfields.u_aux[BUMP_X][grid_ind1] * cos(Omega * fields_ptr->time);
            // a2 = a * a;

	          derivs_ptr->u[U][grid_ind1] =  - globals.dfields.du[U_X][grid_ind1] * a + (- 0.5 * a_x + a_t / a) * (u+v)/2.;
            derivs_ptr->u[V][grid_ind1] =  + globals.dfields.du[V_X][grid_ind1] * a + (  0.5 * a_x + a_t / a) * (u+v)/2.;

            #ifdef DISSIPATION
                // printf("%d \n",globals.dfields.du[DISS_U][grid_ind1]);
                // printf("%d \n",globals.dfields.du[DISS_V][grid_ind1]);
                derivs_ptr->u[U][grid_ind1] = derivs_ptr->u[U][grid_ind1] + sigma*globals.dfields.du[DISS_U][grid_ind1];
                derivs_ptr->u[V][grid_ind1] = derivs_ptr->u[V][grid_ind1] + sigma*globals.dfields.du[DISS_V][grid_ind1];
            #endif // DISSIPATION
        }
    }

    #ifndef PERIODIC

        /*----------------------------- BOUNDARY CONDITIONS --------------------------*/

        /* Boundary conditions at x=xi and x=xf */

        /************************* x = xi ******************************/

        {
            int grid_ind1 = ni_1;

            nx = nx_10;

            dU = sin(omega * fields_ptr->time ) - fields_ptr->u[U][grid_ind1];
            dV = 0.0;

            derivs_ptr->u[U][grid_ind1] =  derivs_ptr->u[U][grid_ind1] + LL*dU;
            derivs_ptr->u[V][grid_ind1] =  derivs_ptr->u[V][grid_ind1] + LL*dV;
        }

        /************************* x = xf ******************************/

        {
            int grid_ind1 = nf_1-1;

            nx = nx_11;

            dU = 0.0;
            dV = - fields_ptr->u[V][grid_ind1];

            derivs_ptr->u[U][grid_ind1] =  derivs_ptr->u[U][grid_ind1] + LL*dU;
            derivs_ptr->u[V][grid_ind1] =  derivs_ptr->u[V][grid_ind1] + LL*dV;
        }

    #endif // Periodic

    /* printf("In eq\n"); */
    /* exit(0); */
}
/*************************************************************************************/
