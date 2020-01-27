/*********************************************************************
*                                                                    *
* This is the Header file for several derivs functions               *
*                                                                    *
*                                                                    *
*********************************************************************/


#include <math.h>
#include "structs_1d.h"    /* place where the used structures are defined */
  
#ifdef DERIVS_1D_H
#else 
#define DERIVS_1D_H

/**********************************************************************
*                                                                     *
* derivD_1d -- take first derivatives of fields using d-nd order with *
*            no boundary condition, in boundary is second order       *
*                                                                     *
* Parameters:                                                         *
*       field    -- pointer to field to derivate                      *
*     dfield    -- ptr to field where fst. derivative is stored       *
*                                                                     *
* Returns: nothing                                                    *
*                                                                     *
**********************************************************************/

extern void derivD_1d(struct GRID_PAR *grid, 
	     FLOAT *field,
	     FLOAT *dfield);


/**********************************************************************
*                                                                     *
* derivQ_1d -- take first derivatives of fields using 4-rth order with*
*            no boundary condition, in boundary is only second order  *
*                                                                     *
* Parameters:                                                         *
*       field    -- pointer to field to derivate                      *
*     dfield    -- ptr to field where fst. derivative is stored       *
*                                                                     *
* Returns: nothing                                                    *
*                                                                     *
**********************************************************************/


extern void derivQ_1d(struct GRID_PAR *grid, 
		 FLOAT *field,
		 FLOAT *dfield);


/*********************************************************************
*                                                                    *
* derivQ_3_3d -- take first derivatives of fields using 4-rth order  *
*             with  no boundary condition, in boundary is only       *
*             third order                                            *
*                                                                    *
* Parameters:                                                        *
*       field    -- pointer to field to derivate                     *
*     dfield    -- ptr to field where fst. derivative is stored      *
*                                                                    *
* Returns: nothing                                                   *
*                                                                    *
*********************************************************************/


extern void derivQ_3_1d(struct GRID_PAR *grid, 
			FLOAT *field,
			FLOAT *dfield);

/*********************************************************************
*                                                                    *
* derivqq -- take derivatives of fields to fourth order, but not     *
*             hermitian                                              *
*                                                                    *
* Returns: nothing                                                   *
*                                                                    *
*********************************************************************/



extern void derivQQ_1d(struct GRID_PAR *grid, 
		       FLOAT *field,
		       FLOAT *dfield);



/*********************************************************************
*                                                                    *
*  deriv_strand_third_order_boundaries_fourth_interior --            *
*      takes first derivatives of fields with 6-th order accuracy    *
*       in the interior, and third order at and close to boundaries. *
*       This is one of Strand's operator, with diagonal scalar       *
*       product. There are more comments in the function itself      * 
*       in derivs_3d.c                                               *
*                                                                    *
* Parameters:                                                        *
*       field    -- pointer to field to derivate                     *
*     dfield    -- ptr to field where fst. derivative is stored      *
*                                                                    *
* Returns: nothing                                                   *
*                                                                    *
*********************************************************************/


extern void  deriv_strand_third_order_boundaries_sixth_interior_1d(struct GRID_PAR *grid,
					  FLOAT *field,
					  FLOAT *dfield
					  );




/*********************************************************************
*                                                                    *
*  deriv_strand_fourth_order_boundaries_eight_interior --            *
*      takes first derivatives of fields with 8-th order accuracy    *
*       in the interior, and fourth order at and close to boundaries.*
*       This is one of Strand's operator, with diagonal scalar       *
*       product. There are more comments in the function itself      * 
*       in derivs_2d.c                                               *
*                                                                    *
* Parameters:                                                        *
*       field    -- pointer to field to derivate                     *
*     dfield    -- ptr to field where fst. derivative is stored      *
*                                                                    *
* Returns: nothing                                                   *
*                                                                    *
*  WARNING: There are two versions of it: each one of them is choosen*
* Using de definitions:                                              *
* NON_OPTIMIZED_EIGHT_ORDER_OPERATOR                                 *
* OPTIMIZED_EIGHT_ORDER_OPERATOR                                     *
*********************************************************************/


extern void  deriv_strand_fourth_order_boundaries_eight_interior_1d(struct GRID_PAR *grid,
					  FLOAT *field,
					  FLOAT *dfield
					  );
/* ----------------- Periodic derivatives ------------------------- */


extern void derivD_Per_1d(struct GRID_PAR *grid,
			  FLOAT *field,
			  FLOAT *dfield
			  );

extern void derivQ_Per_1d(struct GRID_PAR *grid,
			  FLOAT *field,
			  FLOAT *dfield
			  );

extern void deriv_6_Per_1d(struct GRID_PAR *grid,
			   FLOAT *field,
			   FLOAT *dfield
			   );
 
extern void deriv_8_Per_1d(struct GRID_PAR *grid,
			   FLOAT *field,
			   FLOAT *dfield
			   );

/* ------------------ DISSIPATIVE OPERATORS ----------------------- */

/*********************************************************************
*                                                                    *
* diss_KO_4_D_1d -- Kreiss Oliger dissipation of 4th order and       *
*                   second order precision                           *
*                                                                    *
*********************************************************************/

extern void diss_KO4_D_1d(struct GRID_PAR *grid,  
			   FLOAT *field,
			   FLOAT *dfield);

/**********************************************************************
*                                                                     *
* diss_KO4 --  Dissipative fourth order operator, for any of          *
*               Strand's difference operators with diagonal metrics   *
*               Does not destroy the accuracy, but only for second    *
*               order (and first at boundaries) derivatives it is of  *
*               Kreiss-Oliger type.                                   *
*                                                                     *
*                                                                     *
* Parameters:                                                         *
*       field    -- pointer to field to derivate                      *
*     dfield    -- ptr to field where fst. derivative is stored       *
*                                                                     *
* Returns: nothing                                                    *
*                                                                     *
**********************************************************************/

void diss_KO4_1d(struct GRID_PAR *grid,
		 FLOAT *field,
		 FLOAT *dfield);

/**********************************************************************
*                                                                     *
* diss_KO6 --  Dissipative sixth order operator, for any of           *
*               Strand's difference operators with diagonal metrics   *
*               Does not destroy the accuracy, but only for fourth    *
*               order (and second at boundaries) derivatives it is of *
*               Kreiss-Oliger type.                                   *
*                                                                     *
* Parameters:                                                         *
*       field    -- pointer to field to derivate                      *
*     dfield    -- ptr to field where fst. derivative is stored       *
*                                                                     *
* Returns: nothing                                                    *
*                                                                     *
**********************************************************************/

void diss_KO6_1d(struct GRID_PAR *grid,
		 FLOAT *field,
		 FLOAT *dfield);

/**********************************************************************
*                                                                     *
* diss_KO8 --  Dissipative eight order operator, for any of           *
*               Strand's difference operators with diagonal metrics   *
*               Does not destroy the accuracy, but only for sixth     *
*               order (and third at boundaries) derivatives it is of  *
*               Kreiss-Oliger type.                                   *
*                                                                     *
* Parameters:                                                         *
*       field    -- pointer to field to derivate                      *
*     dfield    -- ptr to field where fst. derivative is stored       *
*                                                                     *
* Returns: nothing                                                    *
*                                                                     *
**********************************************************************/

void diss_KO8_1d(struct GRID_PAR *grid,
		 FLOAT *field,
		 FLOAT *dfield);

/******************* PERIODIC VERSIONS *******************************/

void diss_KO4_Per_1d(struct GRID_PAR *grid,
		 FLOAT *field,
		 FLOAT *dfield);


void diss_KO6_Per_1d(struct GRID_PAR *grid,
		 FLOAT *field,
		 FLOAT *dfield);



void diss_KO8_Per_1d(struct GRID_PAR *grid,
		 FLOAT *field,
		 FLOAT *dfield);

#endif
