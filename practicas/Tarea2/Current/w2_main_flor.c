#include "first_macro_1d.h"  /* Where global parameters are defined */
#include "structs_1d.h"      /* Where structures are defined */
#include "derivs_1d.h"       /* Where derivatives functions are defined */
#include "gen_1d.h"

/***********************   Global variables   ****************************/


/* ---------------------------  rkc -------------------------------------*/

#ifdef IMAX
// union fields dv_tempF1;        /* derivative temporal value */
 union fields dv_tempF2;     
 union fields dv_tempF3; 
 union fields dv_tempF4;  
 union fields dv_tempS1; 
 union fields dv_tempS2;  
 union fields dv_tempS3;  
// union fields dv_tempS4;  
 union fields v_temp0; 
 union fields v_temp1;  /* intermediate value of v    */
 
// union fields *dv_tempF1_ptr = &dv_tempF1;        /* derivative temporal value */
 union fields *dv_tempF2_ptr = &dv_tempF2;     
 union fields *dv_tempF3_ptr = &dv_tempF3; 
 union fields *dv_tempF4_ptr = &dv_tempF4;  
 union fields *dv_tempS1_ptr = &dv_tempS1; 
 union fields *dv_tempS2_ptr = &dv_tempS2;  
 union fields *dv_tempS3_ptr = &dv_tempS3;
// union fields *dv_tempS4_ptr = &dv_tempS4; 
 union fields *v_temp0_ptr = &v_temp0; 
 union fields *v_temp1_ptr = &v_temp1;  /* intermediate value of v    */

#else // IMAX

union fields dv_temp;        /* derivative temporal value         */
union fields dv_sum;         /* sum values of deriv               */
union fields v_temp;         /* intermediate value of v           */

union fields *dv_temp_ptr = &dv_temp;    /* derivative temporal value        */
union fields *dv_sum_ptr = &dv_sum;      /* sum          values of deriv     */
union fields *v_temp_ptr = &v_temp;      /* intermediate value of v          */

#endif // IMAX


/* ---------------------------- wave -----------------------------------*/

/* temporal space first derivative values */

struct field_array_derivs dfields;   
struct field_array_aux aux_fields;




/********************************************************************/



int main() {

  /* variable declarations */

  struct GRID_PAR grd;
  struct GRID_PAR *grd_ptr = &grd;
  FLOAT h;
  union fields y; 


  /*   Ploting names */

  struct PLOT_PAR winname; 
  struct PLOT_PAR *winname_ptr = &winname; 


  /*   Initial data parameters */

  struct INI_PAR init_parameters;
  struct INI_PAR *init_par_ptr = &init_parameters;
  
  
  /*   Function parameters  */


  struct FUNCTION_PAR equation_parameters;
  struct FUNCTION_PAR *equation_par_ptr = &equation_parameters;

  /* Parameters coming from first_macro */


#ifdef FILE_INPUT
  FILE *file_data_ptr;
#endif







  /* Get data from web page or data file */
 



INPUT_FUNCTION(grd_ptr, equation_par_ptr, 
		 init_par_ptr, winname_ptr);

 printf("out of input function\n");

file_data_ptr = winname.input_data_file_ptr;
PRINT_MACRO_VALUE(RKX)
PRINT_MACRO_VALUE(DERIV)
PRINT_MACRO_VALUE(FF)
PRINT_MACRO_VALUE(ADISCO)
#ifdef DISSIPATION
PRINT_MACRO_VALUE(DISS)
#endif

/* Array of names for plots */

 sprintf(winname.name[0],"PHI");
 sprintf(winname.name[1],"PHI_T");

  sprintf(winname.window_name[0], "%s_PHI_%d"
	  ,winname.output_file_name,grd.n_gridpts_1);
  sprintf(winname.window_name[1], "%s_PSI_%d"
	  ,winname.output_file_name,grd.n_gridpts_1);


  winname.initial_x= grd.initial_x;
  winname.final_x = grd.final_x;


  winname.n_plots = N_PLOTS;

/* Relation between fields and plot names */

  winname.pointers[0] = PHI;  
  winname.pointers[1] = PHI_T;

/* Number of gridpoints in plots */

  winname.grid_plot_pts_1 = N_GRID_PLOT_PTS_1;  
    
  /* Factor between gridpts and gridplotpts */
 
  winname.factor_1 = FACTOR_1; 
  /* Initial/Final value for time coordinate */

  winname.initial_time = grd.initial_time;       
  winname.final_time = grd.final_time;   

  

  /* Open output file (some times used only for compatibility) */
printf("opening file\n");

	winname_ptr = ADISCO('O',  &winname, &grd, &y); 
//  ADISCO_POINT('O', &winname, &grd, &y);

  /*     creates initial data                            */

inidat(&y,grd_ptr,&init_parameters);

/* write initial data to file                          */ 

/* plot data */
      
	
 
		//#ifdef SDF

	winname_ptr = ADISCO('P',  &winname, &grd, &y); 

		//#endif
	
//    ADISCO_POINT('P',  &winname, &grd, &y); 


  /*     creates potential function */

		// exit(0); 

/* inipot(pot_ptr,&pot_parameters, &grd, winname_ptr); */




  /* makes initial time-interval */

 
  h = grd.time_step; // to be used in the equations.

  /* sends input data to file / browser */

fprintf(OUT_PUT,"<li> Total number of Time Steps = %f </br>\n",(double)(grd.data_steps*grd.int_steps));

fprintf(OUT_PUT,"<li> Number of Time Steps per unit time = %f </br>\n",1.0/h);

fprintf(OUT_PUT,"<li> Time_Step / Space_Step_x = (h/(xf-xi)*(n_gridpts-1)) = %f </br>\n",h*(double)grd.n_gridpts_1/(grd.final_x-grd.initial_x));
fprintf(OUT_PUT,"</ul>%c</br>",10);

fflush(stdout);

/* send input data to the screen */

#ifdef FILE_INPUT
printf("Total number of Time Steps = %f \n",(double)(grd.data_steps*grd.int_steps));

printf("Number of Time Steps per unit time = %f \n",1.0/h);

printf("Time_Step / Space_Step_x= (h/(xf-xi)*n_gridpts) = %f \n",h*(double)grd.n_gridpts_1/(grd.final_x-grd.initial_x));
printf("\n");

fflush(stdout);
#endif

  	 norm_Energy(&grd, &y);

  /* Take data_steps */

  {long int k_outer;
  for (k_outer=1; k_outer<= grd.data_steps; k_outer++) {
/*       printf("h = %f\n",h);  */
/*       printf("time = %f\n",y.a.time); */

#ifdef IMAX
integ(&y,grd_ptr,equation_par_ptr,FF,FS,FI,RKX); 
#else
integ(&y,grd_ptr,equation_par_ptr,FF,RKX); 
#endif //IMAX

/* 	printf("time after integ in main = %f",y.a.time); */
/* printf("Out of integ \n");  */


//      printf("...");
      fflush(stdout);
      /* Do pointwise output */
//      ADISCO_POINT('P',  &winname, &grd, &y); 
      /* Do 1d output */
      if ((k_outer%grd.factor_1d_steps)==0){
			  //#ifdef SDF
		  winname_ptr = ADISCO('P', &winname, &grd, &y);  
			  //#endif
		  
		  norm_Energy(&grd, &y);
      }

	 


/* printf("§</br>\n");  */ 
//printf("r"); 
/* printf("%c",7); */ 
fflush(stdout);
  }



  }
	
//	 winname_ptr = ADISCO('P', &winname, &grd, &y); // printing the last value


fprintf(OUT_PUT,"<ul>%c</br>",10);
fprintf(OUT_PUT,"<li> Execution time = %u secs. ", (unsigned)(clock()/CLOCKS_PER_SEC));
fprintf(OUT_PUT,"</ul>%c</br>",10);

#ifdef FILE_INPUT
printf("\n");
printf("Execution time = %u secs. ", (unsigned)(clock()/CLOCKS_PER_SEC));
printf("\n");
#endif

/* close output file */
  winname_ptr = ADISCO('C',  &winname, &grd, &y); 

#ifdef FILE_INPUT
  fclose(winname.input_data_file_ptr);
#endif
printf("%c",7);
printf("finishing \n");
return(0);
}


void norm_Energy(struct GRID_PAR *grid_1d_ptr,
		union fields *fields_ptr)
		{
		    int ni_1 = (*grid_1d_ptr).start_grid_1; 
			int nf_1 = (*grid_1d_ptr).final_grid_1; 

			
			FLOAT xi = (*grid_1d_ptr).initial_x;
			FLOAT xf = (*grid_1d_ptr).final_x;
			
			FLOAT N=0.0;
			
			int i;

DERIV(grid_1d_ptr,(struct field *)&(*fields_ptr).a.u[PHI][0],
         (struct field *)&dfields.du[PHI_X][0]);

#ifdef PERIODIC

			for(i=ni_1; i < nf_1; i++){
				N = N + (*fields_ptr).a.u[PHI_T][i]*(*fields_ptr).a.u[PHI_T][i] + dfields.du[PHI_X][i]*dfields.du[PHI_X][i];
			}
#else	

  char macro_value_strg[100];

    GET_MACRO_VALUE(DERIV);

    if (strcmp(macro_value_strg,"derivS_1d")==0) {				
      
      				aux_fields.aux[SIGMA][ni_1] = 7493827.0/25401600.0;
      				aux_fields.aux[SIGMA][ni_1+1] = 5534051.0/3628800.0;
      				aux_fields.aux[SIGMA][ni_1+2] = 104561.0/403200.0;
      				aux_fields.aux[SIGMA][ni_1+3] = 260503.0/145152.0;
      				aux_fields.aux[SIGMA][ni_1+4] = 43237.0/103680.0;
      				aux_fields.aux[SIGMA][ni_1+5] = 514081.0/403200.0;
      				aux_fields.aux[SIGMA][ni_1+6] = 3356179.0/3628800.0;
      				aux_fields.aux[SIGMA][ni_1+7] = 25631027.0/25401600.0;
      for(i=ni_1+8; i < nf_1-8; i++){
				aux_fields.aux[SIGMA][i] = 1.0;
			}
 
       				aux_fields.aux[SIGMA][nf_1-1] = 7493827.0/25401600.0;
      				aux_fields.aux[SIGMA][nf_1-2] = 5534051.0/3628800.0;
      				aux_fields.aux[SIGMA][nf_1-3] = 104561.0/403200.0;
      				aux_fields.aux[SIGMA][nf_1-4] = 260503.0/145152.0;
      				aux_fields.aux[SIGMA][nf_1-5] = 43237.0/103680.0;
      				aux_fields.aux[SIGMA][nf_1-6] = 514081.0/403200.0;
      				aux_fields.aux[SIGMA][nf_1-7] = 3356179.0/3628800.0;
      				aux_fields.aux[SIGMA][nf_1-8] = 25631027.0/25401600.0;
 
    }
    
    else if (strcmp(macro_value_strg,"derivQ_1d")==0) {
      				aux_fields.aux[SIGMA][ni_1] = 4567.0/14400.0;
      				aux_fields.aux[SIGMA][ni_1+1] = 799.0/576.0;
      				aux_fields.aux[SIGMA][ni_1+2] = 913.0/1440.0;
      				aux_fields.aux[SIGMA][ni_1+3] = 1769.0/1440.0;
      				aux_fields.aux[SIGMA][ni_1+4] = 2659.0/2880.0;
      				aux_fields.aux[SIGMA][ni_1+5] = 14543.0/14400.0;


	  for(i=ni_1+6; i < nf_1-6; i++){
				aux_fields.aux[SIGMA][i] = 1.0;
			}
			
			      	aux_fields.aux[SIGMA][nf_1-1] = 4567.0/14400.0;
      				aux_fields.aux[SIGMA][nf_1-2] = 799.0/576.0;
      				aux_fields.aux[SIGMA][nf_1-3] = 913.0/1440.0;
      				aux_fields.aux[SIGMA][nf_1-4] = 1769.0/1440.0;
      				aux_fields.aux[SIGMA][nf_1-5] = 2659.0/2880.0;
      				aux_fields.aux[SIGMA][nf_1-6] = 14543.0/14400.0;
    }
    
    else if (strcmp(macro_value_strg,"derivD_1d")==0) {
	   
	   		aux_fields.aux[SIGMA][ni_1] = 0.5;
	   for(i=ni_1+1; i < nf_1-1; i++){
				aux_fields.aux[SIGMA][i] = 1.0;
			}
			aux_fields.aux[SIGMA][nf_1-1] = 0.5;
    }
    
        else if (strcmp(macro_value_strg,"deriv_strand_fourth_order_boundaries_eight_interior_1d")==0) {
	   
	   			    aux_fields.aux[SIGMA][ni_1] = 1498139.0/5080320.0;
      				aux_fields.aux[SIGMA][ni_1+1] = 1107307.0/725760.0;
      				aux_fields.aux[SIGMA][ni_1+2] = 20761.0/80640.0;
      				aux_fields.aux[SIGMA][ni_1+3] = 1304999.0/725760.0;
      				aux_fields.aux[SIGMA][ni_1+4] = 299527.0/725760.0;
      				aux_fields.aux[SIGMA][ni_1+5] = 103097.0/80640.0;
      				aux_fields.aux[SIGMA][ni_1+6] = 670091.0/725760.0;
      				aux_fields.aux[SIGMA][ni_1+7] = 5127739.0/5080320.0;      				      				
	   
	   for(i=ni_1+8; i < nf_1-8; i++){
				aux_fields.aux[SIGMA][i] = 1.0;
		}
		
	   			    aux_fields.aux[SIGMA][nf_1-1] = 1498139.0/5080320.0;
      				aux_fields.aux[SIGMA][nf_1-2] = 1107307.0/725760.0;
      				aux_fields.aux[SIGMA][nf_1-3] = 20761.0/80640.0;
      				aux_fields.aux[SIGMA][nf_1-4] = 1304999.0/725760.0;
      				aux_fields.aux[SIGMA][nf_1-5] = 299527.0/725760.0;
      				aux_fields.aux[SIGMA][nf_1-6] = 103097.0/80640.0;
      				aux_fields.aux[SIGMA][nf_1-7] = 670091.0/725760.0;
      				aux_fields.aux[SIGMA][nf_1-8] = 5127739.0/5080320.0;
      									
    }

 			
    else {
      printf("check sigma value for DERIV!!!%s!!!",macro_value_strg);
    }		

#ifdef SECOND
			
			for(i=ni_1; i < nf_1; i++){
				N = N + aux_fields.aux[SIGMA][i]*((*fields_ptr).a.u[PHI_T][i]*(*fields_ptr).a.u[PHI_T][i] + dfields.du[PHI_X][i]*dfields.du[PHI_X][i]);	
			}
			
#else //SECOND

			for(i=ni_1; i < nf_1; i++){
				N = N + aux_fields.aux[SIGMA][i]*((*fields_ptr).a.u[PHI_T][i]*(*fields_ptr).a.u[PHI_T][i] + (*fields_ptr).a.u[PHI_X][i]*(*fields_ptr).a.u[PHI_X][i]);	
			}

#endif //SECOND			
			  
#endif
			

			
		printf("Time = %f, norm_Energy = %f \n",(*fields_ptr).a.time, N*(xf-xi)/(double)(nf_1 - ni_1));
			
		}
