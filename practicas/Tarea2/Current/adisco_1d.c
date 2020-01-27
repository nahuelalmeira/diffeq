/*********************************************************************
*                                                                    *
* adisco -- sends data to files                                      *
*                                                                    *
* Parameters:                                                        *
*       inst   -- char with instructions ("OPEN", "PUT", or "CLOSE") *
*       fields -- pointer to field_array with data                   *
*                                                                    *
* Returns: pointer to file array                                     *
*                                                                    *
* There are several versions according to the ploting program:       *
* adisco_graph    --  for the grapher program                        *
* adisco_maple3   --  for maple V v.3                                *
* adisco_maple5   --  for maple V v.5                                *
* adisco_mtv      --  for plotmtv                                    *
* adisco_sv       --  for scivis                                     *
*                                                                    *
*********************************************************************/
#include "first_macro_1d.h"  /* Where global parameters are defined */
#include "structs_1d.h"      /* Where structures are defined */
#include "derivs_1d.h"       /* Where derivatives functions are defined */
#include "adisco_1d.h"
#ifdef SV
#include "java_ser.h"
#include "java_fser.h"
#endif
#ifdef SDF
#include "sdf.h"
#endif
#ifdef PYGRAPH
#include "pygraph.h"
#endif


/* ----------------------------------------------------------------------------------------------------------------- */

/* --------------------------------------------- adisco_sdf ------------------------------------------------------*/

/* ----------------------------------------------------------------------------------------------------------------- */


#ifdef SDF

struct PLOT_PAR *adisco_sdf_1d(char inst,
			     struct PLOT_PAR *plot_ptr,
			     struct GRID_PAR *gri,
			       struct field_array *fields){

  char name[100];         /* used to store part of the plot title. */
  char cnames[50];
  char deltanames[50];
  int ret;
  int dim = 1;
  int dim_vector[1];
  double time;


/* number of gridpoints, data from struct plot */

   int n_grid1 = (*gri).n_grid_pts;



 int npp_1 = (*plot_ptr).grid_plot_pts_1;

/* double time = (double)fields->time; */

/* number of variables, data from struct plot */

const int n_plots = (*plot_ptr).n_plots;

    strcpy(cnames, "x");
    strcpy(deltanames, "/nx/");

dim_vector[0] = npp_1;




 switch (inst) {// switch
 case 'O':{// case 'O'
#ifdef DEBUG
printf("write coordinate values in struct plot (SDF) \n");
#endif

/* form coordinate values for plot */

{
double x_i = (double)(*plot_ptr).initial_x;
double x_f = (double)(*plot_ptr).final_x;

(*plot_ptr).coordinate_values = malloc(((*plot_ptr).grid_plot_pts_1) * sizeof(double));

register int l;
       for (l = 0; l < npp_1; l++){
#ifdef NO_LAST_POINT
			(*plot_ptr).coordinate_values[l] = x_i + (double)l/(double)(npp_1)*(x_f-x_i);}
#else
			(*plot_ptr).coordinate_values[l] = x_i + (double)l/(double)(npp_1-1)*(x_f-x_i);}
#endif
printf("npp_1=%d\n",npp_1);

}
}
   break;

 case 'C':{ // case 'C'
#ifdef DEBUG
printf("nothing done \n");
#endif
 }


   break;

 case 'P':{//  case 'P'


#ifdef DEBUG
   printf("writing data in sdf files \n");
#endif

   time = (double)fields->time;

   {register int field_ind;
   for(field_ind = 0; field_ind < n_plots; field_ind++){

       strcpy(name, "Dump/");

	   {register int p_ind1;
	   for(p_ind1 = 0 ; p_ind1 < npp_1; p_ind1++){
	       (*plot_ptr).plot_field[p_ind1] =
		   (double)fields->u[(*plot_ptr).pointers[field_ind]][p_ind1*(*plot_ptr).factor_1];
	   }}

       strcat(name, &(*plot_ptr).window_name[field_ind][0]);
       /* strcat(name, ".nc"); */

	ret=gft_out_full(name,time,&dim_vector[0],cnames,dim,
                              &((*plot_ptr).coordinate_values[0]),(*plot_ptr).plot_field);
#ifdef DEBUG_ADISCO
printf("file name: %s \n",name);
#endif
		}
	}



#ifdef DEBUG
printf("finishing writing data in sdf files \n");
#endif
 }
   break;
 default:
     {
	 printf("Unknown command in first argument of adisco\n");
	 exit(0);
	 break;
     }
 }
return(plot_ptr);
 }


#else // do nothing

struct PLOT_PAR *adisco_dummy_1d(char inst,
							   struct PLOT_PAR *plot_ptr,
							   struct GRID_PAR *gri,
							   struct field_array *fields){

	return(plot_ptr);
}

#endif



/**********************************************************************
 *                              adisco_point                          *
 **********************************************************************/

struct PLOT_PAR *adisco_txt_point_1d(char inst,
			      struct PLOT_PAR *plot_ptr,
			      struct GRID_PAR *grid,
			      struct field_array *fields)
{/* function adisco_txt_point */


  double time = (double)fields->time;
  int n_gridpts_1 = (*grid).n_grid_pts;
  int n_points = (*plot_ptr).n_points;

/*   static FILE *point_output_file_pointer; */
/*   static char point_output_file_name[40]; */



  /* ---------------------------------------------------------------- */

  switch (inst) { //open, close, or write

  case 'O':{ //open files


    {int point_ind;
    char name[100];
    char name2[100];
    for(point_ind=0; point_ind < n_points; point_ind++){
      strcpy(name, "Dump/");
      strcat(name, (*plot_ptr).output_file_name);
      sprintf(name2,"_u2_%d_%d.dat",point_ind,n_gridpts_1);
      strcat(name,name2);

      if (((*plot_ptr).point_output_file_pointer[point_ind] = fopen(name, "w")) == NULL)
	{
	  fprintf(stderr, "Can't open file for pointwise output.\n");
	  exit(1);
	}
    }}}
    break;

  case 'C': {//close files

    {int point_ind;
    for(point_ind=0; point_ind < n_points; point_ind++){
      fclose((*plot_ptr).point_output_file_pointer[point_ind]);
    }}}
    break;

  case 'P':{//plot, or write to files

    time = (double)fields->time;
    double Ave; /* where the average values is stored */
    int N; /* number of points included in average */
      double x_i = (double)(*plot_ptr).initial_x;
      double x_f = (double)(*plot_ptr).final_x;


    {int point_ind;
    register int gr_ind1;

    for(point_ind=0; point_ind < n_points; point_ind++){

  /* coordinates for point averages */

      double x;
      double x0 = (*plot_ptr).x[point_ind];
      double dx = (*plot_ptr).dx[point_ind];
      double dxp = x0+dx;  double dxm = x0-dx;
      Ave=0.0;
      N=0;


	  for(gr_ind1 = 0 ; gr_ind1 < n_gridpts_1; gr_ind1++){
	    x = x_i + (double)gr_ind1/(double)(n_gridpts_1-1)*(x_f-x_i);
	    if(x <= dxp && x > dxm){
		if(point_ind == 0 || point_ind == 1){
		  Ave = Ave + 2.0 * (double)fields->u[U][gr_ind1] * (double)fields->u[U][gr_ind1];
	    }
	    if(point_ind == 2){
		  Ave = Ave + 2.0 * (double)fields->u[V][gr_ind1] * (double)fields->u[V][gr_ind1];
	    }
		  N=N+1;
	    }}
      fprintf((*plot_ptr).point_output_file_pointer[point_ind], "%8f %10e \n", time, Ave/(double)N);
      fflush((*plot_ptr).point_output_file_pointer[point_ind]);
    }
    }
  }

    break;
  default: {
      printf("Unknown command in first argument of adisco\n");
      exit(0);
  }
      break;
  }
return(plot_ptr);
}





/**********************************************************************
 *                              adisco_aschi                          *
 **********************************************************************/

struct PLOT_PAR *adisco_aschi_1d(char inst,
									 struct PLOT_PAR *plot_ptr,
									 struct GRID_PAR *grid,
									 struct field_array *fields)
{/* function adisco_aschi_1d */


	double time = (double)fields->time;
	int npp_1 = (*plot_ptr).grid_plot_pts_1;
	int n_plots = (*plot_ptr).n_plots;
    char name[100];
    char name2[100];


	/*   static FILE *point_output_file_pointer; */
	/*   static char point_output_file_name[40]; */



	/* ---------------------------------------------------------------- */

	switch (inst) { //open, close, or write

		case 'O':{ //open files

			{register int field_ind;
				for(field_ind = 0; field_ind < n_plots; field_ind++){

					strcpy(name, "Dump/");
					strcat(name, &(*plot_ptr).window_name[field_ind][0]);
					sprintf(name2,".dat");
					strcat(name,name2);

		#ifdef DEBUG_ADISCO
					printf("file name: %s \n",name);
		#endif

				if (((*plot_ptr).output_aschi_file_ptr[field_ind] = fopen(name, "w")) == NULL)
					{
						fprintf(stderr, "Can't open file for pointwise output.\n");
						exit(1);
					}

				}
			}
		}
			break;

		case 'C': {//close files

			{int field_ind;
				for(field_ind = 0; field_ind < n_plots; field_ind++){
					fclose((*plot_ptr).output_aschi_file_ptr[field_ind]);
				}}}
			break;

		case 'P':{//plot, or write to files

			time = (double)fields->time;

			{int field_ind;
				register int p_ind1;


				for(field_ind=0; field_ind < n_plots; field_ind++){
					// TOQUE ESTO PARA USAR GNUPLOT
						//fprintf((*plot_ptr).output_aschi_file_ptr[field_ind], "time = %e\n", time);
						fprintf((*plot_ptr).output_aschi_file_ptr[field_ind], "#@type xy\n");
						fprintf((*plot_ptr).output_aschi_file_ptr[field_ind], "\n");
						//fprintf((*plot_ptr).output_aschi_file_ptr[field_ind], "time = %e\n", time);

						for(p_ind1 = 0 ; p_ind1 < npp_1; p_ind1++){
						fprintf((*plot_ptr).output_aschi_file_ptr[field_ind],"%d,     %d,      %e\n", (*plot_ptr).time_slice, p_ind1,
					    (double)fields->u[(*plot_ptr).pointers[field_ind]][p_ind1*(*plot_ptr).factor_1] );
						}
				}
					//fflush((*plot_ptr).output_aschi_file_ptr[field_ind]);

				}
		}

			break;
		default: {
			printf("Unknown command in first argument of adisco\n");
			exit(0);
		}
			break;
	}
	return(plot_ptr);
}


/**********************************************************************
 *                              adisco_pygraph                         *
 **********************************************************************/

struct PLOT_PAR *adisco_pygraph_1d(char inst,
									 struct PLOT_PAR *plot_ptr,
									 struct GRID_PAR *grid,
									 struct field_array *fields)
{/* function adisco_pygraph_1d */


	float time = (float)fields->time;
	int npp_1 = (*plot_ptr).grid_plot_pts_1;
	int n_plots = (*plot_ptr).n_plots;
	int ts = (*plot_ptr).time_slice;
	int field_ind;
	int ierr;
	float dx = (float)((*plot_ptr).final_x - (*plot_ptr).initial_x) / (float)npp_1;
    char name[100];
    char name2[100];


	/*   static FILE *point_output_file_pointer; */
	/*   static char point_output_file_name[40]; */



	/* ---------------------------------------------------------------- */

	switch (inst) { //WRITE (W) and APPEND (A)

		case 'W':{ //open files

			register int p_ind1;


				for(field_ind=0; field_ind < n_plots; field_ind++){

					strcpy(name, "Dump/");
					strcat(name, &(*plot_ptr).window_name[field_ind][0]);
					sprintf(name2,".pyg");
					strcat(name,name2);



						for(p_ind1 = 0 ; p_ind1 < npp_1; p_ind1++){
						(*plot_ptr).plot_field_pygraph[2*p_ind1] = (float)(dx * p_ind1);
						(*plot_ptr).plot_field_pygraph[2*p_ind1+1] =
					    (float)fields->u[(*plot_ptr).pointers[field_ind]][p_ind1*(*plot_ptr).factor_1];
						}

					ierr = pygwrite(name, 'w', ts, time, npp_1, &(*plot_ptr).plot_field_pygraph[0]);
					printf("ierr = %d\n", ierr);
				}

		}
			break;



		case 'A':{//plot, or write to files

				register int p_ind1;


				for(field_ind=0; field_ind < n_plots; field_ind++){

					strcpy(name, "Dump/");
					strcat(name, &(*plot_ptr).window_name[field_ind][0]);
					sprintf(name2,".pyg");
					strcat(name,name2);



						for(p_ind1 = 0 ; p_ind1 < npp_1; p_ind1++){
						(*plot_ptr).plot_field_pygraph[2*p_ind1] = (float)(dx * p_ind1);
						(*plot_ptr).plot_field_pygraph[2*p_ind1+1] =
					    (float)fields->u[(*plot_ptr).pointers[field_ind]][p_ind1*(*plot_ptr).factor_1];
						}

					ierr = pygwrite(name, 'a', ts, time, npp_1, &(*plot_ptr).plot_field_pygraph[0]);
					printf("ierr = %d\n", ierr);
				}
		}

			break;
		default: {
			printf("Unknown command in first argument of adisco\n");
			exit(0);
		}
			break;
	}
	return(plot_ptr);
}




#ifdef SV

/* ----------------------- adisco_sv ----------------------------------*/



struct PLOT_PAR *adisco_sv_1d(char inst,
		       struct PLOT_PAR *plot_ptr,
		       struct GRID_PAR *gri,
		       struct field_array *fields)

{

/* ADISCOSTR=char *data; */

/* number of gridpoints, data from gen.h */
 int npp_1 = (*plot_ptr).grid_plot_pts_1;

/* number of variables, data from gen.h */
const int n_plots = (*plot_ptr).n_plots;

/* The x axis interval */

  double bb[] = {(*plot_ptr).initial_x, (*plot_ptr).final_x};

 switch (inst) {
   case 'O':{
/*    {int field_ind; */
/*       for(field_ind = 0; field_ind < n_fields; ++field_ind){ */
/*        	sprintf((*plot_ptr).name[field_ind], "field(%d)",field_ind); */
/*       }  */
/*    } */
 }
   break;

 case 'C':{
 }
   break;

 case 'P':{
   {
     register int field_ind;
     for(field_ind = 0; field_ind < n_plots; field_ind++){
       {register int p_ind1;
       for(p_ind1 = 0 ; p_ind1 < npp_1; p_ind1++){
	 (*plot_ptr).plot_field[field_ind][p_ind1] =
	   (double)fields->u[(*plot_ptr).pointers[field_ind]][p_ind1*(*plot_ptr).factor_1];
       }
       }
       java_bbser_(&(*plot_ptr).name[field_ind][0], &fields->time, bb,
		   &npp_1, (*plot_ptr).plot_field[field_ind]);
     }
   }
 }
   break;
 default:
   printf("Unknown command \n");
   break;
 }
return(plot_ptr);
}

/* ----------------------- adisco_fsv ----------------------------------*/



struct PLOT_PAR *adisco_fsv_1d(char inst,
			   struct PLOT_PAR *plot_ptr,
			   struct GRID_PAR *gri,
			   struct field_array *fields)

{

char name[100];

/* ADISCOSTR=char *data; */

/* number of gridpoints, data from gen.h */
 int npp_1 = (*plot_ptr).grid_plot_pts_1;

/* number of variables, data from gen.h */
const int n_plots = (*plot_ptr).n_plots;

/* The x axis interval */

  double bb[] = {(*plot_ptr).initial_x, (*plot_ptr).final_x};

strcpy(name,"Dump/");

 switch (inst) {
   case 'O':{
/*    {int field_ind; */
/*       for(field_ind = 0; field_ind < n_fields; ++field_ind){ */
/*        	sprintf((*plot_ptr).name[field_ind], "field(%d)",field_ind); */
/*       }  */
/*    } */
 }
   break;

 case 'C':{
 }
   break;

 case 'P':{
   {register int field_ind;
     for(field_ind = 0; field_ind < n_plots; field_ind++){
       strcat(name, &(*plot_ptr).name[field_ind][0]);

       {register int p_ind1;
       for(p_ind1 = 0 ; p_ind1 < npp_1; p_ind1++){
	 (*plot_ptr).plot_field[field_ind][p_ind1] =
	   (double)fields->u[(*plot_ptr).pointers[field_ind]][p_ind1*(*plot_ptr).factor_1];
       }
       }

       java_fbbser_(&name[0], &fields->time, bb,
		   &npp_1, (*plot_ptr).plot_field[field_ind]);
       strcpy(name, "Dump/");
     }
   }
 }
   break;

 default:
   printf("Unknown command \n");
   break;
 }
return(plot_ptr);
}

#endif
