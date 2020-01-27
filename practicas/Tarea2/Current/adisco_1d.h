
#ifdef ADISCO_1D_H
#else 
#define ADISCO_1D_H

#include "structs_1d.h"
#include "equation.h"


struct plot_1d {
  int n_plots;            /* Number of plots to be produced */
  char name[N_PLOTS][30]; /* Array of names for plots */
  char window_name[N_PLOTS][100]; /* Array of names for plots */
  //FLOAT *plot_field[N_PLOTS]];  /* Fields to plot */


#ifndef SV
/* changed so for reading with fortran convention */
double *plot_field;
float *plot_field_pygraph;
double *coordinate_values;
#endif
#ifdef SV
    float  plot_field[N_PLOTS][N_GRID_PLOT_PTS_1];
#endif

  int pointers[N_PLOTS];  /* Relation between fields and plot names */
  int grid_plot_pts_1;      /* Number of gridpoints in plots */
  int factor_1;             /* Factor between gridpts and gridplotpts */ 

  FLOAT initial_x;       /* Initial value for space coordinate */
  FLOAT final_x;         /* Final value for space coordinate */

  FLOAT initial_time;       /* Initial value for time coordinate */
  FLOAT final_time;         /* Final value for time coordinate */ 
  int time_slice;
 
  /* variables for ploting specific points */

  int n_points; /* number of points to take */
  FILE *point_output_file_pointer[N_POINTS]; /* files where to write values */
  FLOAT x[N_POINTS];   /* coordinate of points */   
  FLOAT dx[N_POINTS];  /* size of region (+-dx) where points are averaged */



  char JSERHOST_env[100];  /* Host where data is sent when using scivis */
  char input_file_name[100];   /* File name where input data goes */
  char output_file_name[100];   /* File name where output data goes */
  FILE *input_data_file_ptr;        /* pointer to the open file */
  FILE *output_data_file_ptr;
	
		// aschi printing
	FILE *output_aschi_file_ptr[N_PLOTS];
};



struct PLOT_PAR *adisco_sv_1d(char inst, 
				   struct PLOT_PAR *plot_ptr, 
				   struct GRID_PAR *gri, 
				   struct field_array  *fields);

struct PLOT_PAR *adisco_fsv_1d(char inst, 
				   struct PLOT_PAR *plot_ptr, 
				   struct GRID_PAR *gri, 
				   struct field_array  *fields);

struct PLOT_PAR *adisco_sdf_1d(char inst, 
				   struct PLOT_PAR *plot_ptr, 
				   struct GRID_PAR *gri, 
				   struct field_array  *fields);



struct PLOT_PAR *adisco_txt_point_1d(char inst,
			     struct PLOT_PAR *plot_ptr,
			     struct GRID_PAR *gri,
			     struct field_array  *fields);

struct PLOT_PAR *adisco_aschi_1d(char inst,
								 struct PLOT_PAR *plot_ptr,
								 struct GRID_PAR *grid,
								 struct field_array  *fields);
								
struct PLOT_PAR *adisco_pygraph_1d(char inst,
								 struct PLOT_PAR *plot_ptr,
								 struct GRID_PAR *grid,
								 struct field_array  *fields);								 


/*********************************************************************
*                                                                    *
* plot_prep   prepares plotting structures to be sent to files       *
*                                                                    *
*********************************************************************/


void plot_prep(struct PLOT_PAR *plot_ptr, 
		      struct FUNCTION_PAR *function_par, 
		      struct field_array  *y_ptr
		      );

#endif
