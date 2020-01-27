/*********************************************************************
*                                                                    *
* input  -- provides initial values/parameters for different         *
* evolutions                                                         *
* Parameters:                                                        *
* y_a_ptr        -- pointer where to write the initial data          *
* initial_time   -- initial time                                     *
*                                                                    *
* Returns: pointer to field_array where data was writen              *
*                                                                    *
*********************************************************************/
#include "first_macro_1d.h"  /* Where global parameters are defined */
#include "structs_1d.h"      /* Where structures are defined */
#include "input.h"



/*-----------------------------------------------------------------------*/

/*********************************************************************************
 *                                                                               *
 *            function to input values from data file                            *
 *                                                                               *
 *********************************************************************************/
  
void w2_1d_input_file(struct GRID_PAR *grid_par, 
			struct FUNCTION_PAR *function_par, 
			struct INI_PAR *ini_pars,
			struct PLOT_PAR *plot_par)
{
	      
    /* input data structure */
    char file_name[30];
    struct entries entries1;
    FILE *file_data_ptr = (*plot_par).input_data_file_ptr;

    /* Get Web Values */
    parse_filedata(&entries1);

    /* File name */
    strcpy(file_name, get_value_char("file_name", &entries1, plot_par));
    /* strcpy((*plot_par).input_file_name, get_value_char("file_name",&entries1,plot_par)); */
    strcpy((*plot_par).output_file_name,file_name);

    strcpy((*plot_par).input_file_name,"input_data_");
    strcat((*plot_par).input_file_name, file_name  /* strcat((*plot_par).input_file_name,"input_data_"); */ );
    (*plot_par).input_data_file_ptr = fopen((*plot_par).input_file_name,"w");
    file_data_ptr = (*plot_par).input_data_file_ptr;

    fprintf(OUT_PUT,"Content-Type: text/html%c%c",10,10);
    fprintf(OUT_PUT,"<body background=\"http://localhost/~reula/Marble.jpg\">%c",10);
    fprintf(OUT_PUT,"<h1>Programa gr-sf - 1D</h1></br>%c",10);
    fprintf(OUT_PUT,"version 0.1 <p>%c",10); 
    printf("Programa gr-sf - 1D \n");
    printf("version 0.1\n"); 

    fprintf(OUT_PUT,"You submitted the following name/value pairs:<p>%c",10);
    fprintf(OUT_PUT,"<ul>%c",10);

    printf("You submitted the following name/value pairs:\n");
    printf("\n");

    /* grid  values */
    (*grid_par).initial_time = get_value_double("initial_time",&entries1,plot_par);
    (*grid_par).final_time = get_value_double("final_time",&entries1,plot_par);
    (*grid_par).data_steps = get_value_int("data_steps",&entries1,plot_par);
    (*grid_par).int_steps = get_value_int("int_steps",&entries1,plot_par);
    (*grid_par).factor_1d_steps = get_value_int("factor_1d_steps",&entries1,plot_par);
    (*grid_par).time_step = ((*grid_par).final_time - (*grid_par).initial_time);
    (*grid_par).time_step =   (*grid_par).time_step / (FLOAT)((*grid_par).data_steps*(*grid_par).int_steps);
 
    fprintf(OUT_PUT,"<li> <code> time_step = %e</code>%c",(*grid_par).time_step,10);

    (*grid_par).initial_x = get_value_double("initial_x",&entries1,plot_par);
    (*grid_par).final_x = get_value_double("final_x",&entries1,plot_par);

    (*grid_par).n_fields = N_FIELDS;
    (*grid_par).start_grid = get_value_int("start_grid",&entries1,plot_par);
    (*grid_par).final_grid = get_value_int("final_grid",&entries1,plot_par);
    (*grid_par).n_grid_pts = (*grid_par).final_grid - (*grid_par).start_grid;

    printf("grd.n_grid_pts = %d \n",(*grid_par).n_grid_pts);


    fprintf(OUT_PUT,"<li> <code> grd.n_grid_pts = %d</code>%c",(*grid_par).n_grid_pts,10);

  /* Function parameters */
  (*function_par).grid_ptr = grid_par;

  (*function_par).c       = get_value_double("c", &entries1, plot_par);
  (*function_par).s       = get_value_double("s", &entries1, plot_par);
  (*function_par).R10     = get_value_double("R10", &entries1, plot_par);
  (*function_par).R11     = get_value_double("R11", &entries1, plot_par);

  (*function_par).nx_10   = get_value_double("nx_10", &entries1, plot_par);
  (*function_par).nx_11   = get_value_double("nx_11", &entries1, plot_par);

  (*function_par).sigma   = get_value_double("sigma", &entries1, plot_par);
  (*function_par).Omega   = get_value_double("Omega", &entries1, plot_par);
  (*function_par).omega   = get_value_double("omega", &entries1, plot_par);
  (*function_par).lambda  = get_value_double("lambda", &entries1, plot_par);

    /* initial data parameters (used in inidat.c and to compare solutions)*/
        
    (*ini_pars).function_par_ptr = function_par;
    (*ini_pars).grid_ptr = grid_par;

    (*ini_pars).a0  = get_value_double("a0",&entries1,plot_par);
    (*ini_pars).k_a_10  = get_value_double("k_a_10",&entries1,plot_par);
    (*ini_pars).shift_a0  = get_value_double("shift_a0",&entries1,plot_par);

    (*ini_pars).a1  = get_value_double("a1",&entries1,plot_par);
    (*ini_pars).k_a_11  = get_value_double("k_a_11",&entries1,plot_par);
    (*ini_pars).shift_a1  = get_value_double("shift_a1",&entries1,plot_par);

    (*ini_pars).c0  = get_value_double("c0",&entries1,plot_par);
    (*ini_pars).k_c_10  = get_value_double("k_c_10",&entries1,plot_par);
    (*ini_pars).shift_c0  = get_value_double("shift_c0",&entries1,plot_par);

    (*ini_pars).c1  = get_value_double("c1",&entries1,plot_par);
    (*ini_pars).k_c_11  = get_value_double("k_c_11",&entries1,plot_par);
    (*ini_pars).shift_c1  = get_value_double("shift_c1",&entries1,plot_par);

    (*ini_pars).b0  = get_value_double("b0",&entries1,plot_par);
    (*ini_pars).sigma_b0  = get_value_double("sigma_b0",&entries1,plot_par);
    (*ini_pars).c0_1  = get_value_double("c0_1",&entries1,plot_par);

    (*ini_pars).b1  = get_value_double("b1",&entries1,plot_par);
    (*ini_pars).sigma_b1  = get_value_double("sigma_b1",&entries1,plot_par);
    (*ini_pars).c1_1  = get_value_double("c1_1",&entries1,plot_par);

    (*ini_pars).v1  = get_value_double("v1",&entries1,plot_par);

    /* Ploting parameters */
    strcpy((*plot_par).JSERHOST_env, "JSERHOST=");
    strcat((*plot_par).JSERHOST_env, get_value_char("host_name",&entries1,plot_par));
    putenv((*plot_par).JSERHOST_env);

    (*plot_par).factor_1 = get_value_int("plot_factor",&entries1,plot_par);
    (*plot_par).grid_plot_pts_1 = ((*grid_par).n_grid_pts -1) / (*plot_par).factor_1 + 1;

    (*plot_par).n_points = N_POINTS; /* the number of points to plot */

    {
        int point_ind;
        char name[100];
        for(point_ind=0; point_ind < (*plot_par).n_points; point_ind++) { 
            sprintf(name,"x[%d]",point_ind);
            (*plot_par).x[point_ind] = get_value_double(name,&entries1,plot_par);
            sprintf(name,"dx[%d]",point_ind);
            (*plot_par).dx[point_ind] = get_value_double(name,&entries1,plot_par);
        }
    }
}

/****************************************************************************
 * Get values from an entry array and asign them to pointer given           *
 *                                                                          *
 * m -- pointer to int giving the number of entries                         * 
 * name[10] -- name of variable as given in entry                           *
 * entries  -- pointer to entry from where to obtain the pairs name/value   *
 * Return   -- value asigned (in double) to name                            *
 ***************************************************************************/

double get_value_double(char name[MAX_ENTRY], 
			struct entries *entries2,
			struct PLOT_PAR *plot_par
			) {
     double interm_float;
     char entry_name[MAX_ENTRY];
     char entry_val[MAX_ENTRY];
     int j;
     int max = (*entries2).dim+1;
     FILE *file_data_ptr = (*plot_par).input_data_file_ptr;


  for(j=0;j<max;++j){
    strcpy(entry_name,(*entries2).entry[j].name);
    strcpy(entry_val,(*entries2).entry[j].val);
    if(strcmp(name,entry_name) == 0) {
      sscanf(entry_val, "%lf", &interm_float);
      fprintf(OUT_PUT,"<li> <code>%s = %f</code>%c",entry_name,(double)interm_float,10);

printf("%s = %f\n",entry_name,(double)interm_float);

      goto L1;
    }
  }
  fprintf(OUT_PUT,"name %s has no counterpart in the entry given", name);
  exit(0);
 L1: return((double)interm_float);
}

int get_value_int(char name[MAX_ENTRY], 
		  struct entries *entries2,
		  struct PLOT_PAR *plot_par
		  ) {
     int interm_int;
     char entry_name[MAX_ENTRY];
     char entry_val[MAX_ENTRY];
     int j;
     int max = (*entries2).dim+1;
     FILE *file_data_ptr = (*plot_par).input_data_file_ptr;

  for(j=0;j<max;++j){
    strcpy(entry_name,(*entries2).entry[j].name);
    strcpy(entry_val,(*entries2).entry[j].val);
    if(strcmp(name,entry_name) == 0) {
      sscanf(entry_val, "%d", &interm_int);
      fprintf(OUT_PUT,"<li> <code>%s = %d</code>%c",entry_name,(int)interm_int,10);

printf("%s = %d\n",entry_name,(int)interm_int);

      goto L1;
    }
  }
  fprintf(OUT_PUT,"name %s has no counterpart in the entry given", name);
  exit(0);
 L1: return((int)interm_int);
}

char *get_value_char(char name[MAX_ENTRY], 
		     struct entries *entries2,
		     struct PLOT_PAR *plot_par
		     ) 
{
    int j;
    int max = (*entries2).dim+1; 
    FILE *file_data_ptr = (*plot_par).input_data_file_ptr;

    for(j=0; j<max; ++j) {
	      if(strcmp(name,(*entries2).entry[j].name) == 0) {
            printf("%s = %s\n",(*entries2).entry[j].name,(*entries2).entry[j].val);
            goto L1;
        }
    }
    fprintf(OUT_PUT,"name %s has no counterpart in the entry given", name);
    exit(0);
    L1: return((*entries2).entry[j].val);
}

void parse_stain(struct entries *entriesp){
    int m=0;
    register int x;
    int cl;
    cl = atoi(getenv("CONTENT_LENGTH"));
    if (cl > ((MAX_ENTRIES+1) * MAX_ENTRY)+10){
	printf("data names too big or too many of them");
	exit(0);
    }
    for(x=0;cl && (!feof(stdin));x++) {
        m=x;
        (*entriesp).entry[x].val = fmakeword(stdin,'&',&cl);
        plustospace((*entriesp).entry[x].val);
        unescape_url((*entriesp).entry[x].val);
        (*entriesp).entry[x].name = makeword((*entriesp).entry[x].val,'=');
    }
    (*entriesp).dim = m;
  }


/*******************************************************************************
 * This function parses a file called data                                     *
 * Comments on that file uses a # at the begining of the file                  *
 * it asign lines of the form name = value                                     *
 * to entries in the structure entries                                         *
 ******************************************************************************/

/* Lee archivo 'data' y almacena los inputs en el puntero entriesp*/
void parse_filedata(struct entries *entriesp)
{
    FILE *fdata_ptr;
    FILE *fdata_ptr_dummy;

    int m = 0;
    int x;
    int incx = 0;
    int cl = 0;

    fdata_ptr = fopen("data","r");
    fdata_ptr_dummy = fopen("data","r");

    while(getc(fdata_ptr_dummy) != EOF) cl++;
    x = 0;
    while(cl && (!feof(fdata_ptr))) {
    
        m = x;
        (*entriesp).entry[x].val = fmakeword(fdata_ptr,'\n',&cl);
	      incx = ((*entriesp).entry[x].val[0] != '#');

	      if(incx) 
            (*entriesp).entry[x].name = makeword((*entriesp).entry[x].val,'=');
	      x = x + incx;
    }
    (*entriesp).dim = m;
    fclose(fdata_ptr);
    fclose(fdata_ptr_dummy);
}
