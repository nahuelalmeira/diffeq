/*********************************************************************
 *                                                                    *
 * This is the Header file of where general structures are defined    * 
 *                                                                    *
 *                                                                    *
 *********************************************************************/

#ifndef TYPES_3D_H
#define TYPES_3D_H

#include "first_macro_3d.h"     /* place where dimensions and float are defined */

/* ----------------> STRUCTURES <---------------------------------- */

struct field_array {
    /* Value of Field at gridpoint */
    FLOAT *u[N_FIELDS];
    FLOAT time;                 /* time at which the values are taken */
    FLOAT radio;
    FLOAT d_t_radio;
/* #ifdef EVOLVE_R */
/* 	   FLOAT lambda; */
/* 	   FLOAT const_C; */
/* #endif */
};


struct deriv_array {
    FLOAT *du[N_DERIVS];
    FLOAT time;                 /* time at which the values are taken */
};


struct aux_array {
    FLOAT *u_aux[N_AUX];
    FLOAT time;                 /* time at which the values are taken */
};


#ifdef TEMPLATE_R
struct template {
    FLOAT *time;
    FLOAT *radio;
    FLOAT *d_t_radio;
};
#endif


struct s_field_array {
    FLOAT *s_x_0_u[N_FIELDS_MPI];       /* boundary values, these are sides x=xi */
    FLOAT *s_x_1_u[N_FIELDS_MPI];       /* boundary values, these are sides x=xf */

    FLOAT *s_y_0_u[N_FIELDS_MPI];       /* boundary values, these are sides y=yi */
    FLOAT *s_y_1_u[N_FIELDS_MPI];       /* boundary values, these are sides y=yf */

    FLOAT *s_z_0_u[N_FIELDS_MPI];       /* boundary values, these are sides z=zi */
    FLOAT *s_z_1_u[N_FIELDS_MPI];       /* boundary values, these are sides z=zf */

    FLOAT time;                 /* time at which the values are taken */
};


struct gcs {
    int g;
    int c;
    int s;
};


struct grid_3d {
    int grid;                   /* The grid of this processor */
    int n_grids;                /* The total number of grids */
    int n_gridpts_1;            /* Number of gridpoints in coordinate 1 */
    int n_gridpts_2;            /* Number of gridpoints in coordinate 2 */
    int n_gridpts_3;            /* Number of gridpoints in coordinate 3 */
    int start_grid_1;           /* Starting value for grid (ussually = zero) */
    int final_grid_1;           /* Final value for the grid (ussually = n_gridpts) */
    int start_grid_2;           /* Starting value for grid (ussually = zero) */
    int final_grid_2;           /* Final value for the grid (ussually = n_gridpts) */
    int start_grid_3;           /* Starting value for grid (ussually = zero) */
    int final_grid_3;           /* Final value for the grid (ussually = n_gridpts) */
    long int n_grid_pts;             /* total number of gridpoints (n_gridpts_1*n_gridpts_2*n_gridpts_3) */
    int layer;                  /* layer of grids on the z-direction start at 1 */
    int n_layers;               /* total number of layers (count from 1) */

    struct gcs gcs_send[N_GRIDS][DIM][2];       /* here store the adjoint grid to given grid and side */
    struct gcs gcs_recv[N_GRIDS][DIM][2];       /* here store the adjoint grid to given grid and side */

    int zgz;                    /* number of ghost-zones along the z direction */
    /* normals to sides */

    FLOAT nx[DIM][2];
    FLOAT ny[DIM][2];
    FLOAT nz[DIM][2];

    /* if excision is valid. This values are given in first_macro.h */
#ifdef EXCISION
    int exc_grid_1;
    int exc_grid_2;
#endif

    int n_fields;               /* Number of fields     */
    int data_steps;             /* Number of time steps where data is collected */
    int int_steps;              /* Number of internal time steps inbetween data collecting */
    int factor_3d_steps;
    FLOAT initial_time;         /* Initial time */
    FLOAT final_time;           /* Final time */
    FLOAT initial_x;
    FLOAT final_x;
    FLOAT initial_y;
    FLOAT final_y;
    FLOAT initial_z;
    FLOAT final_z;
    FLOAT time_step; // =dt= (final_time - initial_time)/data_steps/int_steps

    int n_rk;        // Runge-Kutta steps
    int apply_gauge; // for applying gauge once in a while. 0 false 1 true.
    int gauge_steps; // number of steps between two gauge applications
};


struct plot_3d {
    struct FUNCTION_PAR *function_ptr;
    int n_plots;                /* Number of plots to be produced */
    char name[N_PLOTS][30];     /* Array of names for plots */
    char window_name[N_PLOTS][200];     /* Array of names for plots */
    /*   FLOAT plot_field[N_PLOTS][N_GRID_PLOT_PTS_1][N_GRID_PLOT_PTS_2][N_GRID_PLOT_PTS_3];  */

    /* Fields to plot */

#ifdef NETCDF
    double *plot_field;
#endif
#ifdef SDF
    /* changed so for reading with fortran convention */
    double *plot_field;
    double *coordinate_values;
#endif
#ifdef SV
    float *plot_field;
#endif

    int pointers[N_PLOTS];      /* Relation between fields and plot names */
    int grid_plot_pts_1;        /* Number of gridpoints in plots */
    int factor_1;               /* Factor between gridpts and gridplotpts */
    int grid_plot_pts_2;        /* Number of gridpoints in plots */
    int factor_2;               /* Factor between gridpts and gridplotpts */
    int grid_plot_pts_3;        /* Number of gridpoints in plots */
    int factor_3;               /* Factor between gridpts and gridplotpts */
    int n_gridpts_plot;

    FLOAT initial_x;            /* Initial value for space coordinate */
    FLOAT final_x;              /* Final value for space coordinate */
    FLOAT initial_y;            /* Initial value for space coordinate */
    FLOAT final_y;              /* Final value for space coordinate */
    FLOAT initial_z;            /* Initial value for space coordinate */
    FLOAT final_z;              /* Final value for space coordinate */
    FLOAT initial_time;         /* Initial value for time coordinate */
    FLOAT final_time;           /* Final value for time coordinate */
#ifdef NETCDF
    int tic;                    /* ploting frame, used to mark time in data */
#endif
#ifdef SV
    char JSERHOST_env[40];      /* Host where data is sent when using scivis */
#endif
    char input_file_name[100];  /* File name where input data goes */
    char output_file_name[100]; /* File name where output data goes */
    FILE *input_data_file_ptr;  /* pointer to the open file */
    FILE *output_data_file_ptr;
  
  //------- CHECKPOINTING --------------
  
  int time_level; // time level
//  int prev_time_level;
  int time_stamp;               // value to use in the name of plot file and also to tell that data should be read from files.
  int prev_time_stamp;

    /* To send finishing mail */

    char mailto[200];
    /* variables for ploting specific points */

    int n_points;               /* number of points to take */
    FILE *point_output_file_pointer[N_POINTS];  /* files where to write values */
    FLOAT x[N_POINTS][DIM];     /* coordinate of points */
    FLOAT dx[N_POINTS][DIM];    /* size of region (+-dx) where points are averaged */
};


struct geo {                    // pointwise metric components and inverses, all the geometry needed
    FLOAT xx;
    FLOAT xy;
    FLOAT xz;
    FLOAT yy;
    FLOAT yz;
    FLOAT zz;

    FLOAT uu_xx;
    FLOAT uu_xy;
    FLOAT uu_xz;
    FLOAT uu_yy;
    FLOAT uu_yz;
    FLOAT uu_zz;

    FLOAT det;

    FLOAT nx;
    FLOAT ny;
    FLOAT nz;
    FLOAT n_u_x;
    FLOAT n_u_y;
    FLOAT n_u_z;
    FLOAT norm_n;
    FLOAT norm2_n;
};


struct entry {
    char *name;                 /* Name of data fields in web        */
    char *val;                  /* Value of data fields given in web */
};


struct entries {
    struct entry entry[MAX_ENTRIES];
    int dim;                    /* Real dimension of entry           */
};


struct INI_PAR {
    struct FUNCTION_PAR *function_par_ptr;
    struct GRID_PAR *grid_ptr;
    /* *********** THE COMMENT BELOW IS OUT OF DATE, UPDATE IT OR REMOVE IT **** */
/* Initial data is given by:                                    */
/*       U[0] = PHI = (a0*sin(k_a_10 x + k_a_20 y + k_a_30 z + shift_a0)   */
/*       + c_0*cos(k_c_10 x + k_c_20 y + k_c_30 z + shift_c0))*            */
/*       b_0*exp(-sigma_b0*((x-c0_1)^2+(y-c0_2)^2+(z-c0_3)^2))             */
/* U[1] = dPHI/dt = (a1*sin(k_a_11 x + k_a_21 y + k_a_31 z + shift_a1)     */
/*       + c_1*cos(k_c_11 x + k_c_21 y + k_c_31 z + shift_c1))*            */
/*       b_1*exp(-sigma_b1*((x-c1_1)^2+(y-c1_2)^2+(z-c1_3)^2))             */
/*       + v1*U[0]_x + v2*U[0]_y + v3*U[0]_z                               */

    int initial_data_type;

    // START-initial_data_type = 1 parameters
    int all_patches;
    int spheric_shell;
    FLOAT r_0;
    FLOAT r_0_out;
    FLOAT r_0_in;
    FLOAT a0;
    FLOAT eex;
    FLOAT eey;
    FLOAT eez;
    FLOAT alpha, beta, lambda, gamma;
    int inidat_gridnumber;

    // END-initial_data_type = 1 parameters

#ifdef EVOLVE_R
    FLOAT initial_radio;
    FLOAT lambda;
    FLOAT const_C;
#endif
    FLOAT k_a_10;
    FLOAT k_a_20;
    FLOAT k_a_30;
    FLOAT shift_a0;

    FLOAT a1;
    FLOAT k_a_11;
    FLOAT k_a_21;
    FLOAT k_a_31;
    FLOAT shift_a1;

    FLOAT c0;
    FLOAT k_c_10;
    FLOAT k_c_20;
    FLOAT k_c_30;
    FLOAT shift_c0;

    FLOAT c1;
    FLOAT k_c_11;
    FLOAT k_c_21;
    FLOAT k_c_31;
    FLOAT shift_c1;

    FLOAT b0;
    FLOAT sigma_b0;
    FLOAT c0_1;
    FLOAT c0_2;
    FLOAT c0_3;

    FLOAT b1;
    FLOAT sigma_b1;
    FLOAT c1_1;
    FLOAT c1_2;
    FLOAT c1_3;

    FLOAT v1;
    FLOAT v2;
    FLOAT v3;


    FLOAT sigma_multipole;

};

struct FUNCTION_PAR {
    struct GRID_PAR *grid_ptr;

    FLOAT mm;                   /* mass of field */
    FLOAT m;                    /* mass of BH */
    FLOAT a;                    /* angular moment of BH */
    FLOAT sigma;                /* Parameter for dissipation */
    FLOAT delta;

    int initial_data_type;      /* form of shift function inside horizon */

    FLOAT R10;                  /* Boundary condition at x=0, m = R0*p */
    FLOAT R11;                  /* Boundary condition at x=1, p = R1*m */
    FLOAT R20;                  /* Boundary condition at y=0, m = R0*p */
    FLOAT R21;                  /* Boundary condition at y=1, p = R1*m */

    /* corners */
    FLOAT R_C_00;               /* Corner condition at x=0,y=0,z=0 */
    FLOAT R_C_01;               /* Corner condition at x=0,y=1,z=0 */
    FLOAT R_C_10;               /* Corner condition at x=1,y=0,z=0 */
    FLOAT R_C_11;               /* Corner condition at x=1,y=1,z=0 */
};

struct FLOAT3 {
    FLOAT x;
    FLOAT y;
    FLOAT z;
};
typedef struct FLOAT3 FLOAT3;

struct uint3 {
    unsigned int x;
    unsigned int y;
    unsigned int z;
};
typedef struct uint3 uint3;

typedef enum { X = 0, Y, Z } axis;

typedef enum { XX = 0, XY, XZ, YY, YZ, ZZ, YX, ZX, ZY } plane;

typedef enum { X_XX = 0, X_XY, X_XZ, X_YY, X_YZ, X_ZZ,
               Y_XX, Y_XY, Y_XZ, Y_YY, Y_YZ, Y_ZZ,
               Z_XX, Z_XY, Z_XZ, Z_YY, Z_YZ, Z_ZZ
             } coord3_6;

typedef enum { XX_X = 0, XX_Y, XX_Z,
               XY_X, XY_Y, XY_Z,
               XZ_X, XZ_Y, XZ_Z,
               YY_X, YY_Y, YY_Z,
               YZ_X, YZ_Y, YZ_Z,
               ZZ_X, ZZ_Y, ZZ_Z,
               YX_X, YX_Y, YX_Z,
               ZX_X, ZX_Y, ZX_Z,
               ZY_X, ZY_Y, ZY_Z
             } coord9_3;

typedef enum { XI = 0, XF, YI, YF, ZI, ZF } boundary;

#endif
