#ifndef __INPUT_H
#define __INPUT_H

#include "first_macro_1d.h"
#include "structs_1d.h"
#include "equation.h"
#include "inidat.h"
#include "adisco_1d.h"


struct entry{
  char *name;    /* Name of data fields in web        */
  char *val;     /* Value of data fields given in web */
};

struct entries{
  struct entry entry[MAX_ENTRIES];
  int dim;      /* Real dimension of entry           */
};



/****************************************************************************
 * INPUT_FUNCTION -- This function inputs data to the program which is      *
 * generated via a wave page, the input data structures are given above     *
 *                                                                          *
 ****************************************************************************/
void INPUT_FUNCTION(struct GRID_PAR *grid_parameters,
                    struct FUNCTION_PAR *function_parameters, 
                    struct INI_PAR *ini_parameters,
                    struct PLOT_PAR *plot_parameters);

/****************************************************************************
 * Get values from an entry array and asign them to pointer given           *
 *                                                                          *
 * m -- pointer to int giving the number of entries                         *  
 * name[10] -- name of variable as given in entry                           *
 * entries  -- pointer to entry from where to obtain the pairs name/value   *
 * Return value of name                                                     *
 ***************************************************************************/

FLOAT get_value_double(char name[MAX_ENTRY], struct entries *entries, struct PLOT_PAR *plot_par);
int get_value_int(char name[MAX_ENTRY], struct entries *entries, struct PLOT_PAR *plot_par);
char *get_value_char(char name[MAX_ENTRY], struct entries *entries, struct PLOT_PAR *plot_par);
void parse_stain(struct entries *entriesp);
void parse_filedata(struct entries *entriesp);



/* This are utilities functions to handle data given via web                *
 *                                                                          *
 ****************************************************************************/
 
extern  char *makeword(char *line, char stop);
extern  char *fmakeword(FILE *f, char stop, int *len);
extern  char x2c(char *what);
extern  void unescape_url(char *url);
extern  void plustospace(char *str);


/*************************************************************************/
#endif
