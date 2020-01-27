#define _POSIX_C_SOURCE 200112L
#include <stdlib.h>
#include <sys/time.h>
#include <fenv.h>


#include "first_macro_1d.h"  /* Where global parameters are defined */
#include "structs_1d.h"      /* Where structures are defined */
#include "derivs_1d.h"       /* Where derivatives functions are defined */
//#include "gen_1d.h"
#include "globals.h"
#include "inidat.h"
#include "input.h"
#include "equation.h"
#include "integ.h"
#include "rkc.h"
#include "adisco_1d.h"

// 4M hugepage boundary
#define HUGEPAGE_SIZE (1 << 22)

/***********************   Global variables   ****************************/
struct globals globals;

/***********************   Helper functions   ******************************/

/* Allocate a set of fields at once in huge pages if possible */
void alloc_field_data(size_t nfields, size_t field_elems, FLOAT ** data) {
    FLOAT * memblock;

    size_t alloc_size = nfields * field_elems * sizeof(FLOAT);
    // Request a hugepage-aligned chunk of memory
    if (posix_memalign((void **) &memblock, HUGEPAGE_SIZE, alloc_size) != 0) {
        fprintf(stderr, "out of memory in posix_memaling\n");
        exit(1);
    }

    // Copy pointers to each field
    for (size_t i = 0; i < nfields; ++i) {
        data[i] = &memblock[i * field_elems];
    }
}

/* Free a set of fields allocated with alloc_field_data */
void free_field_data(FLOAT ** data) {
    free(*data);
}

void * safe_malloc(size_t size) {
    void * rv = malloc(size);
    if (rv == NULL) {
        fprintf(stderr, "out of memory in safe_malloc\n");
        exit(1);
    }
    return rv;
}

/***********************   Global variables   ****************************/

void norm_L2(struct GRID_PAR *grid_1d_ptr,
		struct field_array  *fields_ptr);

void norm_Energy(struct GRID_PAR *grid_1d_ptr,
		struct field_array  *fields_ptr);

int main() {

	/* variable declarations */

	struct GRID_PAR grd;
	struct GRID_PAR *grd_ptr = &grd;
	FLOAT h;
	struct field_array y;


	/*   Ploting names */
	struct PLOT_PAR plot;
	struct PLOT_PAR *plot_ptr = &plot;


	/*   Initial data parameters */
	struct INI_PAR init_parameters;
	struct INI_PAR *init_par_ptr = &init_parameters;


	/*   Function parameters  */
	struct FUNCTION_PAR equation_parameters;
	struct FUNCTION_PAR *equation_par_ptr = &equation_parameters;

	/* Parameters coming from first_macro */
	FILE *file_data_ptr;

	/* Get data from web page or data file */

	INPUT_FUNCTION(grd_ptr, equation_par_ptr,
                   init_par_ptr, plot_ptr);

	printf("out of input function\n");

	file_data_ptr = plot_ptr->input_data_file_ptr;
	
	PRINT_MACRO_VALUE(RKX)
	PRINT_MACRO_VALUE(DERIV)
	PRINT_MACRO_VALUE(FF)
	PRINT_MACRO_VALUE(ADISCO)
	
	#ifdef DISSIPATION
		PRINT_MACRO_VALUE(DISS)
	#endif


    /* ------------------------------------ Allocate memory ------------------------------------ */

    /* Allocation #1:------->  auxiliary fields for rkc */
    alloc_field_data(N_FIELDS, grd.n_grid_pts, globals.dv_temp.u);
    alloc_field_data(N_FIELDS, grd.n_grid_pts, globals.dv_sum.u);
    alloc_field_data(N_FIELDS, grd.n_grid_pts, globals.v_temp.u);

    /* Allocation #2: ------> The fields */
    alloc_field_data(N_FIELDS, grd.n_grid_pts, y.u);

    /* Allocation #3:-------> their derivatives (used in evaluating the function) */
    alloc_field_data(N_DERIVS, grd.n_grid_pts, globals.dfields.du);

    /* Allocation #4:--------->  other auxiliary fields */
    alloc_field_data(N_AUX, grd.n_grid_pts, globals.auxfields.u_aux);


	#if defined (SDF) || defined (ASCHI)
    	plot_ptr->plot_field = safe_malloc((plot_ptr->grid_plot_pts_1) * sizeof(double));
	#endif
	#ifdef PYGRAPH
		plot_ptr->plot_field_pygraph = safe_malloc((plot_ptr->grid_plot_pts_1) * sizeof(float)*2);
	#endif


    /* Array of names for plots */

    sprintf(plot_ptr->name[0],"U");
    sprintf(plot_ptr->name[1],"V");

    sprintf(plot_ptr->window_name[0], "%s_U_%d",
            plot_ptr->output_file_name, grd.n_grid_pts);
    sprintf(plot_ptr->window_name[1], "%s_V_%d",
	        plot_ptr->output_file_name, grd.n_grid_pts);

    plot_ptr->initial_x = grd.initial_x;
    plot_ptr->final_x   = grd.final_x;

    plot_ptr->n_plots = N_PLOTS;

    /* Relation between fields and plot names */
    plot_ptr->pointers[0] = U;
    plot_ptr->pointers[1] = V;

    /* Number of gridpoints in plots */

    //plot_ptr->grid_plot_pts_1 = N_GRID_PLOT_PTS_1;

    /* Factor between gridpts and gridplotpts */

    //plot_ptr->factor_1 = FACTOR_1;
    /* Initial/Final value for time coordinate */

    plot_ptr->initial_time = grd.initial_time;
    plot_ptr->final_time   = grd.final_time;

    /* Open output file (some times used only for compatibility) */
    printf("opening file\n");

    #ifdef SDF
		plot_ptr = ADISCO('O',  &plot, &grd, &y);
    #endif

    #ifdef ASCHI
		plot_ptr = adisco_aschi_1d('O',  &plot, &grd, &y);
    #endif
    //  ADISCO_POINT('O', &plot, &grd, &y);

    /* creates initial data */
    inidat(&y, grd_ptr, &init_parameters);

    /* write initial data to file                          */

    /* plot data */
	plot_ptr->time_slice = 0;

    #ifdef SDF
    	plot_ptr = ADISCO('P',  &plot, &grd, &y);
    #endif

    // ADISCO_POINT('P',  &plot, &grd, &y);


    #ifdef ASCHI
		plot_ptr = adisco_aschi_1d('P',  &plot, &grd, &y);
    //	plot_ptr = adisco_aschi_1d('C',  &plot, &grd, &y);
    #endif

    #ifdef PYGRAPH
		plot_ptr = adisco_pygraph_1d('W',  &plot, &grd, &y);
    #endif

    /* creates potential function */



    /* inipot(pot_ptr,&pot_parameters, &grd, plot_ptr); */


    #ifdef BUMP_TEST
        exit(0);
    #endif

    /* makes initial time-interval */

    h = grd.time_step; // to be used in the equations.

    /* sends input data to file / browser */

    fprintf(OUT_PUT, "<li> Total number of Time Steps = %f </br>\n",
	        (double)(grd.data_steps*grd.int_steps));

    fprintf(OUT_PUT, "<li> Number of Time Steps per unit time = %f </br>\n", 1.0/h);

    fprintf(OUT_PUT, 
        "<li> Time_Step / Space_Step_x = (h/(xf-xi)*(n_gridpts-1)) = %f </br>\n",
		h*(double)grd.n_grid_pts/(grd.final_x-grd.initial_x));
    fprintf(OUT_PUT,"</ul>%c</br>", 10);
    fflush(stdout);

    /* send input data to the screen */
	printf("Total number of Time Steps = %f \n",(double)(grd.data_steps*grd.int_steps));
	printf("Number of Time Steps per unit time = %f \n",1.0/h);
	printf("Time_Step / Space_Step_x= (h/(xf-xi)*n_gridpts) = %f \n",h*(double)grd.n_grid_pts/(grd.final_x-grd.initial_x));
	printf("\n");
	fflush(stdout);

  	norm_Energy(&grd, &y);

    /* Take data_steps */

    {
		long int k_outer;
        for (k_outer=1; k_outer<= grd.data_steps; k_outer++) {

            integ(&y, grd_ptr, equation_par_ptr, FF, RKX);
            
			fflush(stdout);

            if ((k_outer%grd.factor_1d_steps)==0){
			    #ifdef SDF
			  	    plot_ptr = ADISCO('P', &plot, &grd, &y);
			    #endif
			    plot_ptr->time_slice = k_outer;
                #ifdef ASCHI
			        plot_ptr = adisco_aschi_1d('P', &plot, &grd, &y);
				#endif
				#ifdef PYGRAPH
					plot_ptr = adisco_pygraph_1d('A',  &plot, &grd, &y);
				#endif
		        norm_Energy(&grd, &y);
            }
			fflush(stdout);
  		}
  	}

    #ifdef SDF
        plot_ptr = ADISCO('P', &plot, &grd, &y);
    #endif

	norm_Energy(&grd, &y);

	fprintf(OUT_PUT,"<ul>%c</br>",10);
	fprintf(OUT_PUT,"<li> Execution time = %u secs. ", (unsigned)(clock()/CLOCKS_PER_SEC));
	fprintf(OUT_PUT,"</ul>%c</br>",10);

	printf("\n");
	printf("Execution time = %u secs. ", (unsigned)(clock()/CLOCKS_PER_SEC));
	printf("\n");

	/* close output file */
	#ifdef SDF
		plot_ptr = ADISCO('C',  &plot, &grd, &y);
	#endif

	#ifdef ASCHI
		plot_ptr = adisco_aschi_1d('C',  &plot, &grd, &y);
	#endif

    fclose(plot_ptr->input_data_file_ptr);

    printf("%c",7);
    printf("finishing \n");
    return(0);
}

/* ------------------------------------------------------------------------------------------------*/
/* ------------------------------------------------------------------------------------------------*/
/* ------------------------------------------------------------------------------------------------*/
/* ------------------------------------------------------------------------------------------------*/
/* ------------------------------------------------------------------------------------------------*/
/* ------------------------------------------------------------------------------------------------*/
/* ------------------------------------------------------------------------------------------------*/
/* ------------------------------------------------------------------------------------------------*/
/* ------------------------------------------------------------------------------------------------*/

void norm_L2(struct GRID_PAR *grid_1d_ptr,
		struct field_array *fields_ptr)
{
	int ni_1 = (*grid_1d_ptr).start_grid;
	int nf_1 = (*grid_1d_ptr).final_grid;

	FLOAT xi = (*grid_1d_ptr).initial_x;
	FLOAT xf = (*grid_1d_ptr).final_x;

	FLOAT N = 0.0;

	int i;

	#ifdef PERIODIC
		for(i=ni_1; i<nf_1; i++){
			N = N + fields_ptr->u[U][i]*fields_ptr->u[U][i];
		}
	#else
		for(i=ni_1+1; i < nf_1-1; i++){
			N = N + fields_ptr->u[U][i]*fields_ptr->u[U][i];
		}
		i = ni_1;
		N = N + 0.5*(fields_ptr->u[U][i]*fields_ptr->u[U][i]);
		i = nf_1-1;
		N = N + 0.5*(fields_ptr->u[U][i]*fields_ptr->u[U][i]);

		nf_1 = nf_1-1; // to make the integral correct
	#endif

	printf("Time = %f, L2 norm = %f \n",fields_ptr->time, N*(xf-xi)/(double)(nf_1 - ni_1));
}

void norm_Energy(struct GRID_PAR *grid_1d_ptr,
		struct field_array *fields_ptr)
{
	int ni_1 = (*grid_1d_ptr).start_grid;
	int nf_1 = (*grid_1d_ptr).final_grid;

	FLOAT xi = (*grid_1d_ptr).initial_x;
	FLOAT xf = (*grid_1d_ptr).final_x;

	FLOAT N=0.0;

	int i;

	#ifdef PERIODIC
		for(i=ni_1; i < nf_1; i++){
			N = N + fields_ptr->u[U][i]*fields_ptr->u[U][i];
		}
	#else

  		char macro_value_strg[100];

      	for(i=ni_1; i < nf_1; i++){
			globals.auxfields.u_aux[SIGMA][i] = 1.0;
		}

    	GET_MACRO_VALUE(DERIV);

    	if (strcmp(macro_value_strg,"derivS_1d")==0) {

			globals.auxfields.u_aux[SIGMA][ni_1] = 7493827.0/25401600.0;
			globals.auxfields.u_aux[SIGMA][ni_1+1] = 5534051.0/3628800.0;
			globals.auxfields.u_aux[SIGMA][ni_1+2] = 104561.0/403200.0;
			globals.auxfields.u_aux[SIGMA][ni_1+3] = 260503.0/145152.0;
			globals.auxfields.u_aux[SIGMA][ni_1+4] = 43237.0/103680.0;
			globals.auxfields.u_aux[SIGMA][ni_1+5] = 514081.0/403200.0;
			globals.auxfields.u_aux[SIGMA][ni_1+6] = 3356179.0/3628800.0;
			globals.auxfields.u_aux[SIGMA][ni_1+7] = 25631027.0/25401600.0;

			globals.auxfields.u_aux[SIGMA][nf_1-1] = 7493827.0/25401600.0;
			globals.auxfields.u_aux[SIGMA][nf_1-2] = 5534051.0/3628800.0;
			globals.auxfields.u_aux[SIGMA][nf_1-3] = 104561.0/403200.0;
			globals.auxfields.u_aux[SIGMA][nf_1-4] = 260503.0/145152.0;
			globals.auxfields.u_aux[SIGMA][nf_1-5] = 43237.0/103680.0;
			globals.auxfields.u_aux[SIGMA][nf_1-6] = 514081.0/403200.0;
			globals.auxfields.u_aux[SIGMA][nf_1-7] = 3356179.0/3628800.0;
			globals.auxfields.u_aux[SIGMA][nf_1-8] = 25631027.0/25401600.0;
	    } 
		else if (strcmp(macro_value_strg,"derivQ_1d")==0) {
			globals.auxfields.u_aux[SIGMA][ni_1] = 4567.0/14400.0;
			globals.auxfields.u_aux[SIGMA][ni_1+1] = 799.0/576.0;
			globals.auxfields.u_aux[SIGMA][ni_1+2] = 913.0/1440.0;
			globals.auxfields.u_aux[SIGMA][ni_1+3] = 1769.0/1440.0;
			globals.auxfields.u_aux[SIGMA][ni_1+4] = 2659.0/2880.0;
			globals.auxfields.u_aux[SIGMA][ni_1+5] = 14543.0/14400.0;

			globals.auxfields.u_aux[SIGMA][nf_1-1] = 4567.0/14400.0;
			globals.auxfields.u_aux[SIGMA][nf_1-2] = 799.0/576.0;
			globals.auxfields.u_aux[SIGMA][nf_1-3] = 913.0/1440.0;
			globals.auxfields.u_aux[SIGMA][nf_1-4] = 1769.0/1440.0;
			globals.auxfields.u_aux[SIGMA][nf_1-5] = 2659.0/2880.0;
			globals.auxfields.u_aux[SIGMA][nf_1-6] = 14543.0/14400.0;
    	} 
		else if (strcmp(macro_value_strg,"derivD_1d")==0) {

	   		globals.auxfields.u_aux[SIGMA][ni_1] = 0.5;

			globals.auxfields.u_aux[SIGMA][nf_1-1] = 0.5;
    	}
        else if (strcmp(macro_value_strg, 
		         "deriv_strand_fourth_order_boundaries_eight_interior_1d")==0) {

			globals.auxfields.u_aux[SIGMA][ni_1] = 1498139.0/5080320.0;
			globals.auxfields.u_aux[SIGMA][ni_1+1] = 1107307.0/725760.0;
			globals.auxfields.u_aux[SIGMA][ni_1+2] = 20761.0/80640.0;
			globals.auxfields.u_aux[SIGMA][ni_1+3] = 1304999.0/725760.0;
			globals.auxfields.u_aux[SIGMA][ni_1+4] = 299527.0/725760.0;
			globals.auxfields.u_aux[SIGMA][ni_1+5] = 103097.0/80640.0;
			globals.auxfields.u_aux[SIGMA][ni_1+6] = 670091.0/725760.0;
			globals.auxfields.u_aux[SIGMA][ni_1+7] = 5127739.0/5080320.0;

			globals.auxfields.u_aux[SIGMA][nf_1-1] = 1498139.0/5080320.0;
			globals.auxfields.u_aux[SIGMA][nf_1-2] = 1107307.0/725760.0;
			globals.auxfields.u_aux[SIGMA][nf_1-3] = 20761.0/80640.0;
			globals.auxfields.u_aux[SIGMA][nf_1-4] = 1304999.0/725760.0;
			globals.auxfields.u_aux[SIGMA][nf_1-5] = 299527.0/725760.0;
			globals.auxfields.u_aux[SIGMA][nf_1-6] = 103097.0/80640.0;
			globals.auxfields.u_aux[SIGMA][nf_1-7] = 670091.0/725760.0;
			globals.auxfields.u_aux[SIGMA][nf_1-8] = 5127739.0/5080320.0;
	    }
    	else {
      		printf("check sigma value for DERIV!!!%s!!!", macro_value_strg);
    	}

		for(i=ni_1; i < nf_1; i++){
			N = N + globals.auxfields.u_aux[SIGMA][i]*fields_ptr->u[U][i]*fields_ptr->u[U][i];
		}
	#endif

	printf("Time = %f, norm_Energy = %1.12e \n",fields_ptr->time, N*(xf-xi)/(double)(nf_1 - ni_1));
}
