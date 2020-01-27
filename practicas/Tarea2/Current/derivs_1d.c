/*********************************************************************
*                                                                    *
* derivs_1d.c -- take derivatives of fields using fourth order method  *
*                                                                    *
* Parameters:                                                        *
*        grid    -- pointer to grid values structure                 *
*       field    -- pointer to field to derivate                     *
*      dfield    -- pointer to field where derivative is to be stored*
*  coordinate    -- pointer to coordinate with respect to which      *
*                   the derivative is taken                          *
* Returns: nothing                                                   *
* Status: tested OK                                                  *
*                                                                    *
*********************************************************************/



#include <math.h>
#include "first_macro_1d.h"
#include "structs_1d.h"    /* place where the used structures are defined */

#ifdef DERIVS_1D_H
#else
#define DERIVS_1D_H


void derivS_1d(struct GRID_PAR *grid,
	     FLOAT *field,
	     FLOAT *dfield)

{/* at function derivD_1d */

    int n1i = (*grid).start_grid;
    int n1f = (*grid).final_grid;

    FLOAT xi = (*grid).initial_x;
    FLOAT xf = (*grid).final_x;

    FLOAT h_1 = (FLOAT)(n1f-n1i-1)/(xf-xi);

    FLOAT h2_1 = 0.5*h_1;


    /* Value at 0 */

    dfield[n1i] = (- field[n1i]
			     + field[n1i+1]
			     )*h_1;


    /* Intermediate Values 1---(n-2) */

    {register int j;
    #pragma omp parallel for
    for (j = n1i+1; j < n1f-1; ++j){
      dfield[j] = (- field[j-1]
			     + field[j+1]
			     )*h2_1;
    }
    }

    /* Value at n-1 */

    dfield[n1f-1] =  (+ field[n1f-1]
				- field[n1f-2]
				)*h_1;
}



void derivQ_1d(struct GRID_PAR *grid,
	     FLOAT *field,
	     FLOAT *dfield)

{/* function derivQ_1d */


  int n1i = (*grid).start_grid;
  int n1f = (*grid).final_grid;

  FLOAT xi = (*grid).initial_x;
  FLOAT xf = (*grid).final_x;

  FLOAT d00 = -1.4117647059;
  FLOAT d01 =  1.7352941176;
  FLOAT d02 = -0.23529411765;
  FLOAT d03 = -0.088235294118;

  FLOAT d10 = -0.5;
  FLOAT d12 = -d10;

  FLOAT d20 = 0.093023255814;
  FLOAT d21 = -0.68604651163;
  FLOAT d23 = -d21;
  FLOAT d24 = -d20;

  FLOAT d30 = 0.030612244898;
  FLOAT d32 = -0.60204081633;
  FLOAT d34 = 0.65306122449;
  FLOAT d35 = -0.081632653061;

  FLOAT d1 = 0.6666666666666;
  FLOAT d2 = -0.0833333333333;

  FLOAT h_1 = (FLOAT)(n1f-n1i-1)/(xf-xi);

  /* Value at 0 */

  dfield[n1i] = (d00*field[n1i]
			   + d01*field[n1i+1]
			   + d02*field[n1i+2]
			   + d03*field[n1i+3]
			   )*h_1;

  /* Value at 1 */

  dfield[n1i+1] = (d10*field[n1i+0]
			     + d12*field[n1i+2]
			     )*h_1;

  /* Value at 2 */

  dfield[n1i+2] = (d20*field[n1i+0]
			     + d21*field[n1i+1]
			     + d23*field[n1i+3]
			     + d24*field[n1i+4]
			     )*h_1;

  /* Value at 3 */

  dfield[n1i+3] = (d30*field[n1i+0]
			     + d32*field[n1i+2]
			     + d34*field[n1i+4]
			     + d35*field[n1i+5]
			     )*h_1;

  /* Intermediate Values 4---(n-5) */

  {register int j;
  #pragma omp parallel for
  for (j = n1i+4; j < n1f-4; ++j){
    dfield[j] = (-d2*field[j-2]
			   - d1*field[j-1]
			   + d1*field[j+1]
			   + d2*field[j+2]
					      )*h_1;
  }
  }

  /* Value at n-1 */

  dfield[n1f-1] = -(d00*field[n1f-1]
			      + d01*field[n1f-2]
			      + d02*field[n1f-3]
			      + d03*field[n1f-4]
			      )*h_1;

  /* Value at n-2 */

  dfield[n1f-2] = -(d10*field[n1f-1]
			      + d12*field[n1f-3]
			      )*h_1;

  /* Value at n-3 */

  dfield[n1f-3] = -( d20*field[n1f-1]
			       + d21*field[n1f-2]
			       + d23*field[n1f-4]
			       + d24*field[n1f-5]
			       )*h_1;

  /* Value at n-4 */

  dfield[n1f-4] = -(d30*field[n1f-1]
			      + d32*field[n1f-3]
			      + d34*field[n1f-5]
			      + d35*field[n1f-6]
			      )*h_1;
}






void derivQ_3_1d(struct GRID_PAR *grid,
	     FLOAT *field,
	     FLOAT *dfield)

{/* function derivQ_3_1d */


    int n1i = (*grid).start_grid;
    int n1f = (*grid).final_grid;

    FLOAT xi = (*grid).initial_x;
    FLOAT xf = (*grid).final_x;

    FLOAT d00 = -1.833333333333333333333333;
    FLOAT d01 =  3.0;
    FLOAT d02 =  -1.5;
    FLOAT d03 =  0.333333333333333333333333;

    FLOAT d10 = -0.38942207148531184298;
    FLOAT d11 = -0.26953763903486946056;
    FLOAT d12 =  0.63903793765926293838;
    FLOAT d13 =  0.094332736084546377480;
    FLOAT d14 = -0.080518371580844513359;
    FLOAT d15 =  0.0061074083572165009295;

    FLOAT d20 =  0.11124996667625322721;
    FLOAT d21 = -0.78615310943278550936;
    FLOAT d22 =  0.19877943763527643222;
    FLOAT d23 =  0.50808067692835148792;
    FLOAT d24 = -0.024137062412656370601;
    FLOAT d25 = -0.0078199093944392672116;

    FLOAT d30 =  0.019051206094885019047822;
    FLOAT d31 =  0.026931104200732614181666;
    FLOAT d32 = -0.633860292039252305642283;
    FLOAT d33 =  0.051772670918649366462688;
    FLOAT d34 =  0.592764606048964306931634;
    FLOAT d35 = -0.054368814269840675877468;
    FLOAT d36 = -0.002290480954138351040607;

    FLOAT d40 = -0.002498706495423627386248;
    FLOAT d41 =  0.005463924453044550084942;
    FLOAT d42 =  0.087024805619019315445041;
    FLOAT d43 = -0.686097670431383548237962;
    FLOAT d44 =  0.018985530480943661987934;
    FLOAT d45 =  0.659895344563505072850627;
    FLOAT d46 = -0.082773228189705424744336;


    FLOAT d1 = 0.6666666666666666666666666;
    FLOAT d2 = -0.0833333333333333333333333;

    FLOAT h_1 = (FLOAT)(n1f-n1i-1)/(xf-xi);


    /* Value at 0 */

    dfield[n1i] = (d00*field[n1i]
			     + d01*field[n1i+1]
			     + d02*field[n1i+2]
			     + d03*field[n1i+3]
			     )*h_1;

    /* Value at 1 */

    dfield[n1i+1] = (d10*field[n1i+0]
			       + d11*field[n1i+1]
			       + d12*field[n1i+2]
			       + d13*field[n1i+3]
			       + d14*field[n1i+4]
			       + d15*field[n1i+5]
			       )*h_1;

    /* Value at 2 */

    dfield[n1i+2] = (d20*field[n1i+0]
			       + d21*field[n1i+1]
			       + d22*field[n1i+2]
			       + d23*field[n1i+3]
			       + d24*field[n1i+4]
			       + d25*field[n1i+5]
			       )*h_1;

    /* Value at 3 */

    dfield[n1i+3] = (d30*field[n1i+0]
			       + d31*field[n1i+1]
			       + d32*field[n1i+2]
			       + d33*field[n1i+3]
			       + d34*field[n1i+4]
			       + d35*field[n1i+5]
			       + d36*field[n1i+6]
			       )*h_1;

    /* Value at 4 */

    dfield[n1i+4] = (d40*field[n1i+0]
			       + d41*field[n1i+1]
			       + d42*field[n1i+2]
			       + d43*field[n1i+3]
			       + d44*field[n1i+4]
			       + d45*field[n1i+5]
			       + d46*field[n1i+6]
			       )*h_1;

    /* Intermediate Values 5---(n-6) */

    {register int j;
    #pragma omp parallel for
    for (j = n1i+5; j < n1f-5; ++j){
      dfield[j] = (-d2*field[j-2]
			     - d1*field[j-1]
			     + d1*field[j+1]
			     + d2*field[j+2]
			     )*h_1;
    }
    }

    /* Value at n-1 */

    dfield[n1f-1] = -(d00*field[n1f-1]
				+ d01*field[n1f-2]
				+ d02*field[n1f-3]
				+ d03*field[n1f-4]
				)*h_1;

    /* Value at n-2 */

    dfield[n1f-2] = -(d10*field[n1f-1]
				+ d11*field[n1f-2]
				+ d12*field[n1f-3]
				+ d13*field[n1f-4]
				+ d14*field[n1f-5]
				+ d15*field[n1f-6]
				)*h_1;

    /* Value at n-3 */

    dfield[n1f-3] = -( d20*field[n1f-1]
				 + d21*field[n1f-2]
				 + d22*field[n1f-3]
				 + d23*field[n1f-4]
				 + d24*field[n1f-5]
				 + d25*field[n1f-6]
				 )*h_1;

    /* Value at n-4 */

    dfield[n1f-4] = -(d30*field[n1f-1]
				+ d31*field[n1f-2]
				+ d32*field[n1f-3]
				+ d33*field[n1f-4]
				+ d34*field[n1f-5]
				+ d35*field[n1f-6]
				+ d36*field[n1f-7]
				)*h_1;
    /* Value at n-5 */

    dfield[n1f-5] = -(d40*field[n1f-1]
				+ d41*field[n1f-2]
				+ d42*field[n1f-3]
				+ d43*field[n1f-4]
				+ d44*field[n1f-5]
				+ d45*field[n1f-6]
					       + d46*field[n1f-7]
				)*h_1;
}







/*********************************************************************
*                                                                    *
* derivqq -- take derivatives of fields to fourth order, but not     *
*             hermitian                                              *
*                                                                    *
* Returns: nothing                                                   *
*                                                                    *
*********************************************************************/



void derivQQ_1d(struct GRID_PAR *grid,
	     FLOAT *field,
	     FLOAT *dfield)

{/* function derivQQ_1d */

    int n1i = (*grid).start_grid;
    int n1f = (*grid).final_grid;

    FLOAT xi = (*grid).initial_x;
    FLOAT xf = (*grid).final_x;




    FLOAT d00 = -2.083333333;
    FLOAT d01 = 4.000000000;
    FLOAT d02 = -3.;
    FLOAT d03 = 1.333333333;
    FLOAT d04 = -.25;

    FLOAT d10 = -.2500000000;
    FLOAT d11 = -.8333333333;
    FLOAT d12 = 1.500000000;
    FLOAT d13 = -.5000000000;
    FLOAT d14 = .08333333333;

    FLOAT dout = .08333333333;
    FLOAT dinn = -.6666666667;

    FLOAT h_1 = (FLOAT)(n1f-n1i-1)/(xf-xi);

    /* Value at 0 */

    dfield[n1i] = (d00*field[n1i+0]
			     + d01*field[n1i+1]
			     + d02*field[n1i+2]
			     + d03*field[n1i+3]
			     + d04*field[n1i+4]
			     )*h_1;

    /* Value at 1 */

    dfield[n1i+1] = (d10*field[n1i+0]
			       + d11*field[n1i+1]
			       + d12*field[n1i+2]
			       + d13*field[n1i+3]
			       + d14*field[n1i+4]
			       )*h_1;



    /* Intermediate Values 3---(n-4) */

    {register int j;
    #pragma omp parallel for
    for (j = n1i+2; j < n1f-2; ++j){
      dfield[j] = (dout*field[j-2]
			     + dinn*field[j-1]
			     - dinn*field[j+1]
			     - dout*field[j+2]
			     )*h_1;
    }
    }

    /* Value at n-1 */

    dfield[n1f-1] = -(d00*field[n1f-1]
				+ d01*field[n1f-2]
				+ d02*field[n1f-3]
				+ d03*field[n1f-4]
				+ d04*field[n1f-5]
				)*h_1;

    /* Value at n-2 */

    dfield[n1f-2] = -(d10*field[n1f-1]
				+ d11*field[n1f-2]
				+ d12*field[n1f-3]
				+ d13*field[n1f-4]
				+ d14*field[n1f-5]
				)*h_1;

}




/* -------------------------------------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------------------------------------- */

void deriv_strand_third_order_boundaries_sixth_interior_1d(struct GRID_PAR *grid,
							 FLOAT *field,
							 FLOAT *dfield)
{


  int n1i = (*grid).start_grid;
  int n1f = (*grid).final_grid;





  FLOAT xi = (*grid).initial_x;
  FLOAT xf = (*grid).final_x;

  /*-----------------------------------------------------------------------------------------------
     I am writing down the closed form expressions (as fractions), for the coefficients of the
     difference operator, the d's (called q's in Strand's paper), so that it is easier to debug, if ever needed.
     The fractions are taken from Strand's paper, I checked the calculations and it seems to me that there
     are no typos there, neither in the difference operator, nor in the scalar product.


     d00=-21600/13649, d01=81763/40947, d02=131/27298, d03=-9143/13649, d04=20539/81894

     d10=-81763/180195, d12=7357/36039, d13=30637/72078, d14=-2328/12013, d15=6611/360390

     d20=-131/54220, d21=-7357/16266, d23=645/2711, d24=11237/32532, d25=-3487/27110

     d30=9143/53590, d31=-30637/64308, d32=-645/5359, d34=-13733/32154, d35=-67/4660, d36=72/5359

     d40=-20539/236310, d41=2328/7877, d42=-11237/47262, d43=-13733/23631, d45=89387/118155,
     d46=-1296/7877, d47=144/7877

     d51=-6611/262806, d52=3487/43801, d53=1541/87602, d54=-89387/131403, d56=32400/43801,
     d57=-6480/43801, d58=720/43801

  ----------------------------------------------------------------------------------------------*/
       //When computed with maple with 32 digits, these numbers are now


  FLOAT d00= -1.5825335189391164187852589933328;
  FLOAT d01= 1.9968007424231323418077026399980;
  FLOAT d02= 0.0047988863653014872884460400029306;
  FLOAT d03= -0.66986592424353432485896402666862;
  FLOAT d04= 0.25079981439421691454807434000049;

  FLOAT d10= -0.45374732928216654180193679069897;
  FLOAT d12= 0.20413995948833208468603457365632;
  FLOAT d13= 0.42505341435666916396126418602070;
  FLOAT d14= -0.19379006076750187297094813951552;
  FLOAT d15= 0.018344016204667166125586170537473;

  FLOAT d20= -0.0024160826263371449649575802286979;
  FLOAT d21= -0.45229312676749047092093938276159;
  FLOAT d23= 0.23791958686831427517521209885651;
  FLOAT d24= 0.34541374646501905815812123447682;
  FLOAT d25= -0.12862412393950571744743637034305;

  FLOAT d30= 0.17061018846799776077626422840082;
  FLOAT d31= -0.47641039995023947253840890713442;
  FLOAT d32= -0.12035827579772345586863220750140;
  FLOAT d34= 0.42710082726876904895191889034024;
  FLOAT d35= -0.014377682403433476394849785407725;
  FLOAT d36= 0.013435342414629595073707781302482;

  FLOAT d40= -0.086915492361728238331005882104016;
  FLOAT d41= 0.29554398882823409927637425415767;
  FLOAT d42= -0.23775972239854428504929964876645;
  FLOAT d43= -0.58114341331302103169565401379544;
  FLOAT d45= 0.75652321103635055647243028225636;
  FLOAT d46= -0.16452964326520248825695061571664;
  FLOAT d47= 0.018281071473911387584105623968516;

  FLOAT d51= -0.025155437851495019139593464380570;
  FLOAT d52= 0.079610054564964270222141047008059;
  FLOAT d53= 0.017590922581676217437958037487729;
  FLOAT d54= -0.68025083141176381056749084876297;
  FLOAT d56= 0.73970913906075203762471176457159;
  FLOAT d57= -0.14794182781215040752494235291432;
  FLOAT d58= 0.016437980868016711947215816990480;

  FLOAT h_1 = (FLOAT)(n1f-n1i-1)/(xf-xi);

  // these are for the intermediate values
  FLOAT d1= 0.750000000000000000000000000000000;
  FLOAT d2= -0.15000000000000000000000000000000;
  FLOAT d3= 0.016666666666666666666666666666667;


#ifdef EXCISION
  printf("************ THIS DIFFERENCE OPERATOR DOES NOT WORK WITH EXCISION !!! **************");
#endif

	  //---------------------->>>> AT AND CLOSE TO LEFT BOUNDARY
	  /* Value at 0 */
	  dfield[n1i] = (d00*field[n1i]
					      + d01*field[n1i+1]
					      + d02*field[n1i+2]
					      + d03*field[n1i+3]
					      + d04*field[n1i+4]
					      )*h_1;
	  /* Value at 1 */
	  dfield[n1i+1] = (d10*field[n1i+0]
						+ d12*field[n1i+2]
						+ d13*field[n1i+3]
						+ d14*field[n1i+4]
						+ d15*field[n1i+5]
						)*h_1;
	  /* Value at 2 */
	  dfield[n1i+2] = (d20*field[n1i+0]
						+ d21*field[n1i+1]
						+ d23*field[n1i+3]
						+ d24*field[n1i+4]
						+ d25*field[n1i+5]
						)*h_1;
	  /* Value at 3 */
	  dfield[n1i+3] = (d30*field[n1i+0]
						+ d31*field[n1i+1]
						+ d32*field[n1i+2]
						+ d34*field[n1i+4]
						+ d35*field[n1i+5]
						+ d36*field[n1i+6]
						)*h_1;
	  /* Value at 4 */
	  dfield[n1i+4] = (d40*field[n1i+0]
						+ d41*field[n1i+1]
						+ d42*field[n1i+2]
						+ d43*field[n1i+3]
						+ d45*field[n1i+5]
						+ d46*field[n1i+6]
						+ d47*field[n1i+7]
						)*h_1;
	  /* Value at 5 */
	  dfield[n1i+5] = (d51*field[n1i+1]
						+ d52*field[n1i+2]
						+ d53*field[n1i+3]
						+ d54*field[n1i+4]
						+ d56*field[n1i+6]
						+ d57*field[n1i+7]
						+ d58*field[n1i+8]
						)*h_1;


	  //------------------>>  INTERMEDIATE VALUES 6-->(n-7)
	  {register int j;
	  #pragma omp parallel for
	  for (j = n1i+6; j < n1f-6; ++j){
	    dfield[j] = (-d3*field[j-3]
					      -d2*field[j-2]
					      - d1*field[j-1]
					      + d1*field[j+1]
					      + d2*field[j+2]
					      + d3*field[j+3]
					      )*h_1;
	  }
	  }


	  //-----------------------> AT AND CLOSE TO RIGHT BOUNDARY
	  /* Value at n-1 */
	  dfield[n1f-1] = -(d00*field[n1f-1]
					      + d01*field[n1f-2]
					      + d02*field[n1f-3]
					      + d03*field[n1f-4]
					      + d04*field[n1f-5]
					      )*h_1;
	  /* Value at n-2 */
	  dfield[n1f-2] = -(d10*field[n1f-1]
						+ d12*field[n1f-3]
						+ d13*field[n1f-4]
						+ d14*field[n1f-5]
						+ d15*field[n1f-6]
						)*h_1;
	  /* Value at n-3 */
	  dfield[n1f-3] = -(d20*field[n1f-1]
						+ d21*field[n1f-2]
						+ d23*field[n1f-4]
						+ d24*field[n1f-5]
						+ d25*field[n1f-6]
						)*h_1;
	  /* Value at n-4 */
	  dfield[n1f-4] = -(d30*field[n1f-1]
						+ d31*field[n1f-2]
						+ d32*field[n1f-3]
						+ d34*field[n1f-5]
						+ d35*field[n1f-6]
						+ d36*field[n1f-7]
						)*h_1;
	  /* Value at n-5 */
	  dfield[n1f-5] = -(d40*field[n1f-1]
						+ d41*field[n1f-2]
						+ d42*field[n1f-3]
						+ d43*field[n1f-4]
						+ d45*field[n1f-6]
						+ d46*field[n1f-7]
						+ d47*field[n1f-8]
						)*h_1;
	  /* Value at n-6 */
	  dfield[n1f-6] = -(d51*field[n1f-2]
						+ d52*field[n1f-3]
						+ d53*field[n1f-4]
						+ d54*field[n1f-5]
						+ d56*field[n1f-7]
						+ d57*field[n1f-8]
						+ d58*field[n1f-9]
						)*h_1;






}



/* ------------------------------------------------------------------------------------------ */


/* ------------------------------------------------------------------------------------------ */






/* ------------------------------------------------------------------------------------------ */


/* ------------------------------------------------------------------------------------------ */

/* ------------------------------------------------------------------------------------------ */


/* ------------------------------------------------------------------------------------------ */


//-----------------------------------------------------------------------------------------

void deriv_strand_fourth_order_boundaries_eight_interior_1d(struct GRID_PAR *grid,
							 FLOAT *field,
							 FLOAT *dfield)
{/* deriv_strand_fourth_order_boundaries_eight_interior */


  int n_gridpts_1 = (*grid).n_grid_pts;

  if (n_gridpts_1<16){
    printf("***************************************************************************************\n");
    printf("error: you are not using enough points for SBP to hold, you need at least 16 in each direction \n");
    printf("and you are using --------> %d  \n", n_gridpts_1);
    printf("***************************************************************************************\n");
    exit(1); }

  int n1i = (*grid).start_grid;
  int n1f = (*grid).final_grid;

#ifdef EXCISION
  printf("************ THIS DIFFERENCE OPERATOR DOES NOT WORK WITH EXCISION !!! **************");
#endif



  FLOAT xi = (*grid).initial_x;
  FLOAT xf = (*grid).final_x;
  /*-----------------------------------------------------------------------------------------------
    The coeffs from Strand's paper for this case are Ok, I checked the calculations and it seems to me that there
     are no typos there, neither in the difference operator, nor in the scalar product.
     ----------------------------------------------------------------------------------------------*/
  //When computed with maple with 32 digits, these numbers are now

#ifdef NON_OPTIMIZED_EIGHT_ORDER_OPERATOR
  FLOAT d00=-1.6955436044318985087498556542484;
  FLOAT d01=2.0610513554928258770826116045752;
  FLOAT d02=0.87789728901434824583477679084963;
  FLOAT d03=-2.5445639556810149125014434575163;
  FLOAT d04=1.6889486445071741229173883954248;
  FLOAT d05=-0.38778972890143482458347767908496;
  FLOAT d06= 0.0;
  FLOAT d07= 0.0;
  FLOAT d08= 0.0;
  FLOAT d09= 0.0;
  FLOAT d010= 0.0;
  FLOAT d011= 0.0;

  FLOAT d10=-0.39835918734972926809517745906661;
  FLOAT d11= 0.0;
  FLOAT d12=-0.44127885642072764523900478066757;
  FLOAT d13=1.8989658393441626095262349706691;
  FLOAT d14=-1.5738365692621829357170143420027;
  FLOAT d15=0.60604617027316423238240764906811;
  FLOAT d16=-0.091537396584686992857446038000302;
  FLOAT d17= 0.0;
  FLOAT d18= 0.0;
  FLOAT d19= 0.0;
  FLOAT d110= 0.0;
  FLOAT d111= 0.0;

  FLOAT d20=-1.0055577090642838920860210778808;
  FLOAT d21=2.6151125596961539353524325346492;
  FLOAT d22= 0.0;
  FLOAT d23=-10.448835244859039473942416960576;
  FLOAT d24=16.854276269429682078398421541815;
  FLOAT d25=-10.479793844227156688020808246231;
  FLOAT d26=2.5403284648744886691035720465902;
  FLOAT d27=-0.075530495849844628805179838367064;
  FLOAT d28= 0.0;
  FLOAT d29= 0.0;
  FLOAT d210= 0.0;
  FLOAT d211= 0.0;

  FLOAT d30=0.41730852995727528198434306408336;
  FLOAT d31=-1.6112948490126556929673253900322;
  FLOAT d32=1.4960581706703734383449080548465;
  FLOAT d33= 0.0;
  FLOAT d34=-1.6738166082885887268879133240715;
  FLOAT d35=1.9668117242490862700533359284822;
  FLOAT d36=-0.64585509260926636725392126737262;
  FLOAT d37=0.050788125033775796726572934064281;
  FLOAT d38= 0.0;
  FLOAT d39= 0.0;
  FLOAT d310= 0.0;
  FLOAT d311= 0.0;

  FLOAT d40=-1.2067978767157806109219504339143;
  FLOAT d41=5.8182409265274916785465083281307;
  FLOAT d42=-10.513925845304986417473772536922;
  FLOAT d43=7.2925946575767793888363987219850;
  FLOAT d44= 0.0;
  FLOAT d45=-3.0385143353798934097204370001146;
  FLOAT d46=2.0057652009112144592418490931814;
  FLOAT d47=-0.34870908370292399111170040001164;
  FLOAT d48=-0.0086536439119010973968957723343805;
  FLOAT d49= 0.0;
  FLOAT d410= 0.0;
  FLOAT d411= 0.0;

  FLOAT d50=0.089446187550400088107295882497877;
  FLOAT d51=-0.72324463225750363106445242685870;
  FLOAT d52=2.1103523865873885758072494834961;
  FLOAT d53=-2.7662054325681783318765966175687;
  FLOAT d54=0.98086385026111691291085453864196;
  FLOAT d55= 0.0;
  FLOAT d56=0.30802293704705996043280348352271;
  FLOAT d57=-0.026238992169135626268229544755535;
  FLOAT d58=0.029797181295285022842565739061272;
  FLOAT d59=-0.0027934857464329708914905380369943;
  FLOAT d510= 0.0;
  FLOAT d511= 0.0;

  FLOAT d60= 0.0;
  FLOAT d61=0.15126303740835199995224529205735;
  FLOAT d62=-0.70834831886017471258878769201994;
  FLOAT d63=1.2577996869081960509841200672744;
  FLOAT d64=-0.89656603854302375846464634405377;
  FLOAT d65=-0.42651843804299217071512177699248;
  FLOAT d66= 0.0;
  FLOAT d67=0.80159349003841766765011518336067;
  FLOAT d68=-0.21661535522787203529072916962025;
  FLOAT d69=0.041260067662451816245853175165761;
  FLOAT d610=-0.0038681313433548577730487351717901;
  FLOAT d611= 0.0;

  FLOAT d70= 0.0;
  FLOAT d71= 0.0;
  FLOAT d72=0.019265719907611002302054245220619;
  FLOAT d73=-0.090478311526126167237971095382715;
  FLOAT d74=0.14258418768974005892265577479665;
  FLOAT d75=0.033235928479719164073418453370319;
  FLOAT d76=-0.73326354624003548802568409455577;
  FLOAT d77= 0.0;
  FLOAT d78=0.79260196355547737511601116983528;
  FLOAT d79=-0.19815049088886934377900279245882;
  FLOAT d710=0.037742950645498922624571960468347;
  FLOAT d711=-0.0035384016230155239960536212939075;
#endif

  // these are the coeffs for the improved operator that we choose not by minimizing the bandwith, but
  // by minimizing the maximum imaginary eigenvalue.

#ifdef OPTIMIZED_EIGHT_ORDER_OPERATOR
  FLOAT d00= -1.6955436044318985087498556542484;
  FLOAT d01= 2.252142137678813514633822362278;
  FLOAT d02= -0.05888015064022764242837280119;
  FLOAT d03= -0.72703849019795003445385686285;
  FLOAT d04= -0.03519446459907925766567721687;
  FLOAT d05= 0.3808994692748803682435341446956;
  FLOAT d06= -0.097708425809176140086689775338;
  FLOAT d07= -0.018676471275362299492904196473;
  FLOAT d08= 0.;
  FLOAT d09= 0.;
  FLOAT d010=0.;
  FLOAT d011= 0.;

  FLOAT d10= -0.4352931378302752275823881078793;
  FLOAT d11= 0.;
  FLOAT d12= 0.099854313212144418846805809049;
  FLOAT d13= 0.48598825799891087114955473054;
  FLOAT d14= -0.04056967339078804101602657017;
  FLOAT d15= -0.1516077655067655130871564977011;
  FLOAT d16= 0.028751917941456163466861493696;
  FLOAT d17= 0.012876087575317328222349142474;
  FLOAT d18= 0.;
  FLOAT d19= 0.;
  FLOAT d110= 0.;
  FLOAT d111= 0.;

  FLOAT d20= 0.06744227386055814358882611857;
  FLOAT d21= -0.59175794358011014241446301559;
  FLOAT d22= 0.;
  FLOAT d23= 0.09922290191545044394136441726;
  FLOAT d24= 1.24445434548753271358155515947;
  FLOAT d25= -1.368619045325369683541255238187;
  FLOAT d26= 0.68543519098309329993738259237;
  FLOAT d27= -0.13617772334115477509341003392;
  FLOAT d28= 0.;
  FLOAT d29= 0.;
  FLOAT d210= 0.;
  FLOAT d211= 0.;

  FLOAT d30= 0.119234324171531019965386581760;
  FLOAT d31= -0.41236675277145806241997120304;
  FLOAT d32= -0.01420667755300961916445912985;
  FLOAT d33= 0.;
  FLOAT d34= -0.117132848377661592077848335517;
  FLOAT d35= 0.6750458046328004849045861337825;
  FLOAT d36= -0.288099428939536862991210465807;
  FLOAT d37= 0.037525578837334631783516418677;
  FLOAT d38= 0.;
  FLOAT d39= 0.;
  FLOAT d310= 0.;
  FLOAT d311= 0.;

  FLOAT d40= 0.025147363295176347088194768037;
  FLOAT d41= 0.14998007970344354042651691946;
  FLOAT d42= -0.77630747812384192410033152270;
  FLOAT d43= 0.51033212364828546341398271274;
  FLOAT d44= 0.;
  FLOAT d45= -0.568060141489748837333529197702;
  FLOAT d46= 0.891763513806768671939424492616;
  FLOAT d47= -0.224201816928182164037362400114;
  FLOAT d48= -0.0086536439119010973968957723343805;
  FLOAT d49= 0.;
  FLOAT d410= 0.;
  FLOAT d411= 0.;

  FLOAT d50= -0.0878569049859194092294958469532;
  FLOAT d51= 0.1809259887937250033140311227933;
  FLOAT d52= 0.275603557814485387547649301144;
  FLOAT d53= -0.949412365700909499468138419805;
  FLOAT d54= 0.1833756882676831850910631088522;
  FLOAT d55= 0.;
  FLOAT d56= 0.42315722086966643064298670184389;
  FLOAT d57= -0.052796880607583149849171168899192;
  FLOAT d58= 0.029797181295285022842565739061272;
  FLOAT d59= -0.0027934857464329708914905380369943;
  FLOAT d510= 0.;
  FLOAT d511= 0.;

  FLOAT d60= 0.031207020141045519688121961545;
  FLOAT d61= -0.047511755865994320174424070761;
  FLOAT d62= -0.191127593117949651614482212113;
  FLOAT d63= 0.561072252375672358928364455971;
  FLOAT d64= -0.398613397284846386535560095569;
  FLOAT d65= -0.58594453589139385546142240382277;
  FLOAT d66= 0.;
  FLOAT d67= 0.81014142855224141198732709437972;
  FLOAT d68= -0.21661535522787203529072916962025;
  FLOAT d69= 0.041260067662451816245853175165761;
  FLOAT d610= -0.0038681313433548577730487351717901;
  FLOAT d611= 0.;

  FLOAT d70= 0.0054565862264050490869367571164;
  FLOAT d71= -0.019463641447689387726897436343;
  FLOAT d72= 0.034735133749982204632490070184;
  FLOAT d73= -0.066851276946818081029475174146;
  FLOAT d74= 0.091674221978406727279476067978;
  FLOAT d75= 0.066875790674993403525413442454852;
  FLOAT d76= -0.74108283592437134573347044379599;
  FLOAT d77= 0.;
  FLOAT d78= 0.79260196355547737511601116983528;
  FLOAT d79= -0.19815049088886934377900279245882;
  FLOAT d710= 0.037742950645498922624571960468347;
  FLOAT d711= -0.0035384016230155239960536212939075;
#endif



  FLOAT h_1 = (FLOAT)(n1f-n1i-1)/(xf-xi);
  // these are for the intermediate values
  FLOAT d1= 0.800000000000000000000000000000000;
  FLOAT d2= -0.2000000000000000000000000000000;
  FLOAT d3= 0.038095238095238095238095238095238;
  FLOAT d4= -0.0035714285714285714285714285714286;



	  //---------------------->>>> AT AND CLOSE TO LEFT BOUNDARY
	  /* Value at 0 */
	  dfield[n1i] = (d00*field[n1i+0]
					      + d01*field[n1i+1]
					      + d02*field[n1i+2]
					      + d03*field[n1i+3]
					      + d04*field[n1i+4]
					      + d05*field[n1i+5]
					      + d06*field[n1i+6]
					      + d07*field[n1i+7]
					      + d08*field[n1i+8]
					      + d09*field[n1i+9]
					      + d010*field[n1i+10]
					      + d011*field[n1i+11]
					      )*h_1;
	  /* Value at 1 */
	  dfield[n1i+1] = (d10*field[n1i+0]
						+ d11*field[n1i+1]
						+ d12*field[n1i+2]
						+ d13*field[n1i+3]
						+ d14*field[n1i+4]
						+ d15*field[n1i+5]
						+ d16*field[n1i+6]
						+ d17*field[n1i+7]
						+ d18*field[n1i+8]
						+ d19*field[n1i+9]
						+ d110*field[n1i+10]
						+ d111*field[n1i+11]
						)*h_1;
	  /* Value at 2 */
	  dfield[n1i+2] = (d20*field[n1i+0]
						+ d21*field[n1i+1]
						+ d22*field[n1i+2]
						+ d23*field[n1i+3]
						+ d24*field[n1i+4]
						+ d25*field[n1i+5]
						+ d26*field[n1i+6]
						+ d27*field[n1i+7]
						+ d28*field[n1i+8]
						+ d29*field[n1i+9]
						+ d210*field[n1i+10]
						+ d211*field[n1i+11]
						)*h_1;
	  /* Value at 3 */
	  dfield[n1i+3] = (d30*field[n1i+0]
						+ d31*field[n1i+1]
						+ d32*field[n1i+2]
						+ d33*field[n1i+3]
						+ d34*field[n1i+4]
						+ d35*field[n1i+5]
						+ d36*field[n1i+6]
						+ d37*field[n1i+7]
						+ d38*field[n1i+8]
						+ d39*field[n1i+9]
						+ d310*field[n1i+10]
						+ d311*field[n1i+11]
						)*h_1;
	  /* Value at 4 */
	  dfield[n1i+4] = (d40*field[n1i+0]
						+ d41*field[n1i+1]
						+ d42*field[n1i+2]
						+ d43*field[n1i+3]
						+ d44*field[n1i+4]
						+ d45*field[n1i+5]
						+ d46*field[n1i+6]
						+ d47*field[n1i+7]
						+ d48*field[n1i+8]
						+ d49*field[n1i+9]
						+ d410*field[n1i+10]
						+ d411*field[n1i+11]
						)*h_1;
	  /* Value at 5 */
	  dfield[n1i+5] = (d50*field[n1i+0]
						+ d51*field[n1i+1]
						+ d52*field[n1i+2]
						+ d53*field[n1i+3]
						+ d54*field[n1i+4]
						+ d55*field[n1i+5]
						+ d56*field[n1i+6]
						+ d57*field[n1i+7]
						+ d58*field[n1i+8]
						+ d59*field[n1i+9]
						+ d510*field[n1i+10]
						+ d511*field[n1i+11]
						)*h_1;

	  /* Value at 6 */
	  dfield[n1i+6] = (d60*field[n1i+0]
						+ d61*field[n1i+1]
						+ d62*field[n1i+2]
						+ d63*field[n1i+3]
						+ d64*field[n1i+4]
						+ d65*field[n1i+5]
						+ d66*field[n1i+6]
						+ d67*field[n1i+7]
						+ d68*field[n1i+8]
						+ d69*field[n1i+9]
						+ d610*field[n1i+10]
						+ d611*field[n1i+11]
						)*h_1;

	  /* Value at 7 */
	  dfield[n1i+7] = (d70*field[n1i+0]
						+ d71*field[n1i+1]
						+ d72*field[n1i+2]
						+ d73*field[n1i+3]
						+ d74*field[n1i+4]
						+ d75*field[n1i+5]
						+ d76*field[n1i+6]
						+ d77*field[n1i+7]
						+ d78*field[n1i+8]
						+ d79*field[n1i+9]
						+ d710*field[n1i+10]
						+ d711*field[n1i+11]
						)*h_1;


	  //------------------>>  INTERMEDIATE VALUES 6-->(n-7)
	  {register int j;
	  #pragma omp parallel for
	  for (j = n1i+8; j < n1f-8; ++j){
	    dfield[j] = (-d4*field[j-4]
					      -d3*field[j-3]
					      -d2*field[j-2]
					      - d1*field[j-1]
					      + d1*field[j+1]
					      + d2*field[j+2]
					      + d3*field[j+3]
					      + d4*field[j+4]
					      )*h_1;
	  }
	  }


	  //-----------------------> AT AND CLOSE TO RIGHT BOUNDARY

	  /* Value at n-1 */
	  dfield[n1f-1] = -(d00*field[n1f-1]
					      + d01*field[n1f-2]
					      + d02*field[n1f-3]
					      + d03*field[n1f-4]
					      + d04*field[n1f-5]
					      + d05*field[n1f-6]
					      + d06*field[n1f-7]
					      + d07*field[n1f-8]
					      + d08*field[n1f-9]
					      + d09*field[n1f-10]
					      + d010*field[n1f-11]
					      + d011*field[n1f-12]
					      )*h_1;
	  /* Value at n-2 */
	  dfield[n1f-2] = -(d10*field[n1f-1]
					      + d11*field[n1f-2]
					      + d12*field[n1f-3]
					      + d13*field[n1f-4]
					      + d14*field[n1f-5]
					      + d15*field[n1f-6]
					      + d16*field[n1f-7]
					      + d17*field[n1f-8]
					      + d18*field[n1f-9]
					      + d19*field[n1f-10]
					      + d110*field[n1f-11]
					      + d111*field[n1f-12]
					      )*h_1;

	  /* Value at n-3 */
	  dfield[n1f-3] = -(d20*field[n1f-1]
					      + d21*field[n1f-2]
					      + d22*field[n1f-3]
					      + d23*field[n1f-4]
					      + d24*field[n1f-5]
					      + d25*field[n1f-6]
					      + d26*field[n1f-7]
					      + d27*field[n1f-8]
					      + d28*field[n1f-9]
					      + d29*field[n1f-10]
					      + d210*field[n1f-11]
					      + d211*field[n1f-12]
					      )*h_1;

	  /* Value at n-4 */
	  dfield[n1f-4] = -(d30*field[n1f-1]
					      + d31*field[n1f-2]
					      + d32*field[n1f-3]
					      + d33*field[n1f-4]
					      + d34*field[n1f-5]
					      + d35*field[n1f-6]
					      + d36*field[n1f-7]
					      + d37*field[n1f-8]
					      + d38*field[n1f-9]
					      + d39*field[n1f-10]
					      + d310*field[n1f-11]
					      + d311*field[n1f-12]
					      )*h_1;

	    /* Value at n-5 */
	    dfield[n1f-5] = -(d40*field[n1f-1]
					      + d41*field[n1f-2]
					      + d42*field[n1f-3]
					      + d43*field[n1f-4]
					      + d44*field[n1f-5]
					      + d45*field[n1f-6]
					      + d46*field[n1f-7]
					      + d47*field[n1f-8]
					      + d48*field[n1f-9]
					      + d49*field[n1f-10]
					      + d410*field[n1f-11]
					      + d411*field[n1f-12]
					      )*h_1;
	  /* Value at n-6 */
	  dfield[n1f-6] = -(d50*field[n1f-1]
					      + d51*field[n1f-2]
					      + d52*field[n1f-3]
					      + d53*field[n1f-4]
					      + d54*field[n1f-5]
					      + d55*field[n1f-6]
					      + d56*field[n1f-7]
					      + d57*field[n1f-8]
					      + d58*field[n1f-9]
					      + d59*field[n1f-10]
					      + d510*field[n1f-11]
					      + d511*field[n1f-12]
					      )*h_1;
	  /* Value at n-7 */
	  dfield[n1f-7] = -(d60*field[n1f-1]
					      + d61*field[n1f-2]
					      + d62*field[n1f-3]
					      + d63*field[n1f-4]
					      + d64*field[n1f-5]
					      + d65*field[n1f-6]
					      + d66*field[n1f-7]
					      + d67*field[n1f-8]
					      + d68*field[n1f-9]
					      + d69*field[n1f-10]
					      + d610*field[n1f-11]
					      + d611*field[n1f-12]
					      )*h_1;

	  /* Value at n-8 */
	  dfield[n1f-8] =  -(d70*field[n1f-1]
					      + d71*field[n1f-2]
					      + d72*field[n1f-3]
					      + d73*field[n1f-4]
					      + d74*field[n1f-5]
					      + d75*field[n1f-6]
					      + d76*field[n1f-7]
					      + d77*field[n1f-8]
					      + d78*field[n1f-9]
					      + d79*field[n1f-10]
					      + d710*field[n1f-11]
					      + d711*field[n1f-12]
					      )*h_1;



}



/* ------------------------------------------------------------------------------------------ */

/* -------------------------------------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------------------------------------- */


/* -----------------------------------------------PERIODIC OPERATORS------------------------------------------ */


void derivD_Per_1d(struct GRID_PAR *grid,
	     FLOAT *field,
	     FLOAT *dfield)

{/* at function derivD_Per_1d */

    int n1i = (*grid).start_grid;
    int n1f = (*grid).final_grid;

    FLOAT xi = (*grid).initial_x;
    FLOAT xf = (*grid).final_x;

    // FLOAT h_1 = (FLOAT)(n1f-n1i-1)/(xf-xi);
    FLOAT h_1 = (FLOAT)(n1f-n1i)/(xf-xi);

    FLOAT h2_1 = 0.5*h_1;


#ifdef EXCISION
  printf("************ THIS DIFFERENCE OPERATOR DOES NOT WORK WITH EXCISION !!! **************");
#endif


    /* Value at 0 */

    dfield[n1i] = (- field[n1f-1]
			     + field[n1i+1]
			     )*h2_1;


    /* Intermediate Values 1---(n-2) */

    {register int j;
    #pragma omp parallel for
    for (j = n1i+1; j < n1f-1; ++j){
      dfield[j] = (- field[j-1]
			     + field[j+1]
			     )*h2_1;
    }
    }

    /* Value at n-1 */

    dfield[n1f-1] =  (+ field[n1i]
				- field[n1f-2]
				)*h2_1;
}



void derivQ_Per_1d(struct GRID_PAR *grid,
	     FLOAT *field,
	     FLOAT *dfield)

{/* function derivQ_Per_1d */

  int n1i = (*grid).start_grid;
  int n1f = (*grid).final_grid;

  FLOAT xi = (*grid).initial_x;
  FLOAT xf = (*grid).final_x;


  FLOAT d1 = 0.6666666666666;
  FLOAT d2 = -0.0833333333333;

  // FLOAT h_1 = (FLOAT)(n1f-n1i-1)/(xf-xi);
  FLOAT h_1 = (FLOAT)(n1f-n1i)/(xf-xi);

#ifdef EXCISION
  printf("************ THIS DIFFERENCE OPERATOR DOES NOT WORK WITH EXCISION !!! **************");
#endif


  /* Value at 0 */

  dfield[n1i] = (- d2*field[n1f-2]
			   - d1*field[n1f-1]
			   + d1*field[n1i+1]
			   + d2*field[n1i+2]
			   )*h_1;

  /* Value at 1 */

  dfield[n1i+1] = (- d2*field[n1f-1]
			   - d1*field[n1i]
			   + d1*field[n1i+2]
			   + d2*field[n1i+3]
			     )*h_1;



  /* Intermediate Values 4---(n-5) */

  {register int j;
  #pragma omp parallel for
  for (j = n1i+2; j < n1f-2; ++j){
    dfield[j] = (-d2*field[j-2]
			   - d1*field[j-1]
			   + d1*field[j+1]
			   + d2*field[j+2]
					      )*h_1;
  }
  }

  /* Value at n-1 */

  dfield[n1f-1] = (-d2*field[n1f-3]
			   - d1*field[n1f-2]
			   + d1*field[n1i]
			   + d2*field[n1i+1]
			      )*h_1;

  /* Value at n-2 */

  dfield[n1f-2] = (-d2*field[n1f-4]
			   - d1*field[n1f-3]
			   + d1*field[n1f-1]
			   + d2*field[n1i]
			      )*h_1;


}



void deriv_6_Per_1d(struct GRID_PAR *grid,
							 FLOAT *field,
							 FLOAT *dfield)
{/* deriv_6_Per_1d */



  int n1i = (*grid).start_grid;
  int n1f = (*grid).final_grid;





  FLOAT xi = (*grid).initial_x;
  FLOAT xf = (*grid).final_x;






  FLOAT h_1 = (FLOAT)(n1f-n1i)/(xf-xi); /* For periodic boundary conditions we take out the 1 */

  // these are for the intermediate values
  FLOAT d1= 0.750000000000000000000000000000000;
  FLOAT d2= -0.15000000000000000000000000000000;
  FLOAT d3= 0.016666666666666666666666666666667;


#ifdef EXCISION
  printf("************ THIS DIFFERENCE OPERATOR DOES NOT WORK WITH EXCISION !!! **************");
#endif

	  //---------------------->>>> AT AND CLOSE TO LEFT BOUNDARY
	  /* Value at 0 */
	  dfield[n1i] = (- d3*field[n1f-3]
				   - d2*field[n1f-2]
				   - d1*field[n1f-1]
				   + d1*field[n1i+1]
				   + d2*field[n1i+2]
				   + d3*field[n1i+3]
				   )*h_1;
	  /* Value at 1 */
	  dfield[n1i+1] = (- d3*field[n1f-2]
				     - d2*field[n1f-1]
				     - d1*field[n1i]
				     + d1*field[n1i+2]
				     + d2*field[n1i+3]
				     + d3*field[n1i+4]
				     )*h_1;
	  /* Value at 2 */
	  dfield[n1i+2] = (- d3*field[n1f-1]
				     - d2*field[n1i]
				     - d1*field[n1i+1]
				     + d1*field[n1i+3]
				     + d2*field[n1i+4]
				     + d3*field[n1i+5]
				     )*h_1;

	  //------------------>>  INTERMEDIATE VALUES 6-->(n-7)
	  {register int j;
	  #pragma omp parallel for
	  for (j = n1i+3; j < n1f-3; ++j){
	    dfield[j] = (- d3*field[j-3]
				   - d2*field[j-2]
				   - d1*field[j-1]
				   + d1*field[j+1]
				   + d2*field[j+2]
				   + d3*field[j+3]
				   )*h_1;
	  }
	  }


	  //-----------------------> AT AND CLOSE TO RIGHT BOUNDARY
	  /* Value at n-1 */
	  dfield[n1f-1] = (- d3*field[n1f-4]
				     - d2*field[n1f-3]
				     - d1*field[n1f-2]
				     + d1*field[n1i]
				     + d2*field[n1i+1]
				     + d3*field[n1i+2]
				     )*h_1;
	  /* Value at n-2 */
	  dfield[n1f-2] = (- d3*field[n1f-5]
				     - d2*field[n1f-4]
				     - d1*field[n1f-3]
				     + d1*field[n1f-1]
				     + d2*field[n1i]
				     + d3*field[n1i+1]
				     )*h_1;
	  /* Value at n-3 */
	  dfield[n1f-3] = (- d3*field[n1f-6]
				     - d2*field[n1f-5]
				     - d1*field[n1f-4]
				     + d1*field[n1f-2]
				     + d2*field[n1f-1]
				     + d3*field[n1i]
				     )*h_1;


}


void deriv_8_Per_1d(struct GRID_PAR *grid,
							 FLOAT *field,
							 FLOAT *dfield)
{/* deriv_8_Per_1d */



  int n1i = (*grid).start_grid;
  int n1f = (*grid).final_grid;





  FLOAT xi = (*grid).initial_x;
  FLOAT xf = (*grid).final_x;






  FLOAT h_1 = (FLOAT)(n1f-n1i)/(xf-xi); /* For periodic boundary conditions we take out the 1 */

  // these are for the intermediate values

  FLOAT d1= 0.800000000000000000000000000000000;
  FLOAT d2= -0.2000000000000000000000000000000;
  FLOAT d3= 0.038095238095238095238095238095238;
  FLOAT d4= -0.0035714285714285714285714285714286;




#ifdef EXCISION
  printf("************ THIS DIFFERENCE OPERATOR DOES NOT WORK WITH EXCISION !!! **************");
#endif

	  //---------------------->>>> AT AND CLOSE TO LEFT BOUNDARY
	  /* Value at 0 */
	  dfield[n1i] = (- d4*field[n1f-4]
				   - d3*field[n1f-3]
				   - d2*field[n1f-2]
				   - d1*field[n1f-1]
				   + d1*field[n1i+1]
				   + d2*field[n1i+2]
				   + d3*field[n1i+3]
				   + d4*field[n1i+4]
				   )*h_1;
	  /* Value at 1 */
	  dfield[n1i+1] = (- d4*field[n1f-3]
				     - d3*field[n1f-2]
				     - d2*field[n1f-1]
				     - d1*field[n1i]
				     + d1*field[n1i+2]
				     + d2*field[n1i+3]
				     + d3*field[n1i+4]
				     + d4*field[n1i+5]
				     )*h_1;
	  /* Value at 2 */
	  dfield[n1i+2] = (- d4*field[n1f-2]
				     - d3*field[n1f-1]
				     - d2*field[n1i]
				     - d1*field[n1i+1]
				     + d1*field[n1i+3]
				     + d2*field[n1i+4]
				     + d3*field[n1i+5]
				     + d4*field[n1i+6]
				     )*h_1;
	  /* Value at 3 */
	  dfield[n1i+3] = (- d4*field[n1f-1]
				     - d3*field[n1i]
				     - d2*field[n1i+1]
				     - d1*field[n1i+2]
				     + d1*field[n1i+4]
				     + d2*field[n1i+5]
				     + d3*field[n1i+6]
				     + d4*field[n1i+7]
				     )*h_1;

	  //------------------>>  INTERMEDIATE VALUES 6-->(n-7)
	  {register int j;
	  #pragma omp parallel for
	  for (j = n1i+4; j < n1f-4; ++j){
	    dfield[j] = (- d4*field[j-4]
				   - d3*field[j-3]
				   - d2*field[j-2]
				   - d1*field[j-1]
				   + d1*field[j+1]
				   + d2*field[j+2]
				   + d3*field[j+3]
				   + d4*field[j+4]
				   )*h_1;
	  }
	  }


	  //-----------------------> AT AND CLOSE TO RIGHT BOUNDARY
	  /* Value at n-1 */
	  dfield[n1f-1] = (- d4*field[n1f-5]
				     - d3*field[n1f-4]
				     - d2*field[n1f-3]
				     - d1*field[n1f-2]
				     + d1*field[n1i]
				     + d2*field[n1i+1]
				     + d3*field[n1i+2]
				     + d4*field[n1i+3]
				     )*h_1;
	  /* Value at n-2 */
	  dfield[n1f-2] = (- d4*field[n1f-6]
				     - d3*field[n1f-5]
				     - d2*field[n1f-4]
				     - d1*field[n1f-3]
				     + d1*field[n1f-1]
				     + d2*field[n1i]
				     + d3*field[n1i+1]
				     + d4*field[n1i+2]
				     )*h_1;
	  /* Value at n-3 */
	  dfield[n1f-3] = (- d4*field[n1f-7]
				     - d3*field[n1f-6]
				     - d2*field[n1f-5]
				     - d1*field[n1f-4]
				     + d1*field[n1f-2]
				     + d2*field[n1f-1]
				     + d3*field[n1i]
				     + d4*field[n1i+1]
				     )*h_1;

	  /* Value at n-4 */
	  dfield[n1f-4] = (- d4*field[n1f-8]
				     - d3*field[n1f-7]
				     - d2*field[n1f-6]
				     - d1*field[n1f-5]
				     + d1*field[n1f-3]
				     + d2*field[n1f-2]
				     + d3*field[n1f-1]
				     + d4*field[n1i]
				     )*h_1;






}







/* -------------------------------------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------------------------------------- */


/* -------------------------------------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------------------------------------- */


/* -----------------------------------------------DISSIPATIVE OPERATORS------------------------------------------ */


/* -------------------------------------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------------------------------------- */


void diss_KO4_D_1d(struct GRID_PAR *grid,
	     FLOAT *field,
	     FLOAT *dfield)

{/* at function diss_KO_4_D_1d */

    int n1i = (*grid).start_grid;
    int n1f = (*grid).final_grid;

    FLOAT xi = (*grid).initial_x;
    FLOAT xf = (*grid).final_x;

    FLOAT h_1 = (FLOAT)(n1f-n1i-1)/(xf-xi);



    /* Value at 0 */

    dfield[n1i] = -2.0*(      field[n1i+0]
			          -2.0* field[n1i+1]
				  +     field[n1i+2]
					)*h_1;

    dfield[n1i+1] = -(-2.0*field[n1i+0]
			        +5.0*field[n1i+1]
				-4.0*field[n1i+2]
				+    field[n1i+3]
					)*h_1;


    /* Intermediate Values 2---(n-3) */

    {register int j;
    #pragma omp parallel for
    for (j = n1i+2; j < n1f-2; ++j){
      dfield[j] = (-     field[j-2]
			     + 4.0*field[j-1]
			     - 6.0*field[j-0]
			     + 4.0*field[j+1]
			     -     field[j+2]
			     )*h_1;
    }
    }


    /* Value at n-1 */



    dfield[n1f-1] = -2.0*(      field[n1f-1]
			            -2.0* field[n1f-2]
				    +     field[n1f-3]
					  )*h_1;

    dfield[n1f-2] = -(-2.0*field[n1f-1]
			        +5.0*field[n1f-2]
				-4.0*field[n1f-3]
				+    field[n1f-4]
					)*h_1;

}




/* ********************************************************************************************************** */
/* ********************************************************************************************************** */
/* ********************************************************************************************************** */
/* ********************************************************************************************************** */

/* NEW DISSIPATION OPERATORS */

/* ********************************************************************************************************** */
/* ********************************************************************************************************** */
/* ********************************************************************************************************** */
/* ********************************************************************************************************** */


//-------------------------- KO-type dissipation, fourth derivative in the interior --------------------------//
// Written by Manuel.
void diss_KO4_1d(struct GRID_PAR *grid,
	      FLOAT *field,
		 FLOAT *dfield)
{// diss_KO4

#ifdef EXCISION
  printf("************ THIS DIFFERENCE OPERATOR DOES NOT WORK WITH EXCISION !!! **************");
#endif


  int n1i = (*grid).start_grid;
  int n1f = (*grid).final_grid;

  FLOAT xi = (*grid).initial_x;
  FLOAT xf = (*grid).final_x;

  FLOAT d2 = 2.0;
  FLOAT d4 = 4.0;
  FLOAT d5 = 5.0;
  FLOAT d6 = 6.0;

  FLOAT p1;


  FLOAT h_1 = (xf-xi)/(FLOAT)(n1f-n1i-1); // this is dx

  int n1r;
  int n1s;


  // in the radial direction
  FLOAT sigma_rm2; // scalar product at gridpoint R Minus 2
  FLOAT sigma_rm1; // scalar product at gridpoint R Minus 1
  FLOAT sigma_sp2; // scalar product at gridpoint S Plus 2
  FLOAT sigma_sp1; // scalar product at gridpoint S Plus 1



  // define what scalar product, interior means (ie which are gridpoints r,s), depending on the derivatives that we are using

  char macro_value_strg[100];
  GET_MACRO_VALUE(DERIV);

  // --------- derivatives on the sphere ------------
  if (strcmp(macro_value_strg,"derivS_1d")==0) { // second order in the interior
    n1r = n1i+2;  n1s = n1f-3;
    p1 = 1.0/h_1;
    sigma_rm2 = 0.5;  sigma_rm1 = 1.0;
    sigma_sp2 = sigma_rm2;  sigma_sp1 = sigma_rm1;
  }
  else if (strcmp(macro_value_strg,"derivQ_1d")==0) {// fourth order in the interior
    n1r = n1i+4;  n1s = n1f-5;
    p1 = h_1;
    sigma_rm2 = 43.0/48.0;  sigma_rm1 = 49.0/48.0;
    sigma_sp2 = sigma_rm2;  sigma_sp1 = sigma_rm1;
  }
  else if (strcmp(macro_value_strg,"deriv_strand_third_order_boundaries_sixth_interior_1d")==0) { // sixth order in the interior
    n1r = n1i+6;  n1s = n1f-7;
    p1 = pow(h_1,3);
    sigma_rm2 = 7877.0/8640.0;  sigma_rm1 = 43801.0/43200.0;
    sigma_sp2 = sigma_rm2;  sigma_sp1 = sigma_rm1;
  }
  else if (strcmp(macro_value_strg,"deriv_strand_fourth_order_boundaries_eight_interior_1d")==0) {// eighth order in the interior
  n1r = n1i+8;  n1s = n1f-9;
  p1 = pow(h_1,5);
  sigma_rm2 = 670091.0/725760.0;  sigma_rm1 = 5127739.0/5080320.0;
  sigma_sp2 = sigma_rm2;  sigma_sp1 = sigma_rm1;
  }
  else {
    printf("I don't know what derivative you are using, will stop here.\n");
    exit(1);
  }


  // initialize to zero
  {register int ind1;
  for (ind1 = n1i; ind1 < n1f-1; ind1++){
	dfield[ind1] =  0.0;
      }
    }


  //------------------------- dissipation in the x-direction ----------------------------//
  // interior
  {register int ind1;
  #pragma omp parallel for
  for (ind1 = n1r; ind1 < n1s+1; ind1++){

	dfield[ind1] =  -(   field[ind1-2]  -
				    d4*field[ind1-1]  +
				    d6*field[ind1+0] -
				    d4*field[ind1+1]  +
				       field[ind1+2])*p1;
      }
    }


  // boundaries
  {register int ind1;

      // left one
      ind1 = n1r-2;
      dfield[ind1] =  -(   field[ind1+0]  -
				  d2*field[ind1+1]  +
				     field[ind1+2])*p1/sigma_rm2;

      ind1 = n1r-1;
      dfield[ind1] =  -( - d2*field[ind1-1]  +
				     d5*field[ind1+0]  -
				     d4*field[ind1+1]  +
				        field[ind1+2])*p1/sigma_rm1;


      // right one
      ind1 = n1s+1;
      dfield[ind1] =  -( - d2*field[ind1+1]  +
				     d5*field[ind1+0]  -
				     d4*field[ind1-1]  +
				        field[ind1-2])*p1/sigma_sp1;

      ind1 = n1s+2;
      dfield[ind1] = -(      field[ind1+0]  -
				    d2*field[ind1-1]  +
				       field[ind1-2])*p1/sigma_sp2;


    }



}


//-------------------------- KO-type dissipation, sixth derivative in the interior --------------------------//
// Written by Manuel.

void diss_KO6_1d(struct GRID_PAR *grid,
	      FLOAT *field,
		 FLOAT *dfield)
{// diss_KO6




  int n1i = (*grid).start_grid;
  int n1f = (*grid).final_grid;

  FLOAT xi = (*grid).initial_x;
  FLOAT xf = (*grid).final_x;



  FLOAT p1;

  FLOAT h_1 = (xf-xi)/(FLOAT)(n1f-n1i-1); // this is dx

  int n1r;
  int n1s;

  // scalar product
  FLOAT sigma_rm3; // scalar product at gridpoint R Minus 3
  FLOAT sigma_rm2; // scalar product at gridpoint R Minus 2
  FLOAT sigma_rm1; // scalar product at gridpoint R Minus 1
  FLOAT sigma_sp3; // scalar product at gridpoint S Plus 3
  FLOAT sigma_sp2; // scalar product at gridpoint S Plus 2
  FLOAT sigma_sp1; // scalar product at gridpoint S Plus 1




  // define what scalar product, interior means (ie which are gridpoints r,s), depending on the derivatives that we are using

  char macro_value_strg[100];
  GET_MACRO_VALUE(DERIV);

  // --------- derivatives on the sphere ------------
  if (strcmp(macro_value_strg,"derivS_1d")==0) { // second order in the interior
    n1r = n1i+3;  n1s = n1f-4;
    p1 = 1.0/pow(h_1,3);
    sigma_rm3 = 0.5; sigma_rm2 = 1.0;  sigma_rm1 = 1.0;
    sigma_sp3 = sigma_rm3; sigma_sp2 = sigma_rm2;  sigma_sp1 = sigma_rm1;
  }

  else if (strcmp(macro_value_strg,"derivQ_1d")==0) {// fourth order in the interior
    n1r = n1i+4;  n1s = n1f-5;
    p1 = 1.0/h_1;
    sigma_rm3 = 59.0/48.0; sigma_rm2 = 43.0/48.0;  sigma_rm1 = 49.0/48.0;
    sigma_sp3 = sigma_rm3; sigma_sp2 = sigma_rm2;  sigma_sp1 = sigma_rm1;
  }

  else if (strcmp(macro_value_strg,"deriv_strand_third_order_boundaries_sixth_interior_1d")==0) { // sixth order in the interior
    n1r = n1i+6;  n1s = n1f-7;
    p1 = h_1;
    sigma_rm3= 5359.0/4320.0; sigma_rm2 = 7877.0/8640.0;  sigma_rm1 = 43801.0/43200.0;
    sigma_sp3 = sigma_rm3; sigma_sp2 = sigma_rm2;  sigma_sp1 = sigma_rm1;
  }

  else if (strcmp(macro_value_strg,"deriv_strand_fourth_order_boundaries_eight_interior_1d")==0) {// eighth order in the interior
  n1r = n1i+8;  n1s = n1f-9;
  p1 = pow(h_1,3);
  sigma_rm3 = 103097.0/80640.0; sigma_rm2 = 670091.0/725760.0;  sigma_rm1 = 5127739.0/5080320.0;
  sigma_sp3 = sigma_rm3; sigma_sp2 = sigma_rm2;  sigma_sp1 = sigma_rm1;
  }

  else {
    printf("I don't know what derivative you are using, will stop here.\n");
    exit(1);
  }

#ifdef EXCISION
  printf("************ THIS DIFFERENCE OPERATOR DOES NOT WORK WITH EXCISION !!! **************");
#endif

  // initialize to zero
  {register int ind1;
  for (ind1 = n1i; ind1 < n1f-1; ind1++){

	dfield[ind1] =  0.0;

      } }


  //------------------------- dissipation in the x-direction ----------------------------//
  // interior
  {register int ind1;
  #pragma omp parallel for
  for (ind1 = n1r; ind1 < n1s+1; ind1++){
	dfield[ind1] =  ( field[ind1+3]
						     - 6.0*field[ind1+2]
						     + 15.0*field[ind1+1]
						     - 20.0*field[ind1+0]
						     + 15.0*field[ind1-1]
						     - 6.0*field[ind1-2]
						     + field[ind1-3])*p1;
      } }

  // boundaries
  {register int ind1;

      // left one
      ind1 = n1r-3;
      dfield[ind1] =  (field[ind1+3]
						  - 3.0*field[ind1+2]
						  + 3.0*field[ind1+1]
						  - field[ind1+0])*p1/sigma_rm3;

      ind1 = n1r-2;
      dfield[ind1] =  (field[ind1+3]
						  - 6.0*field[ind1+2]
						  + 12.0*field[ind1+1]
						  - 10.0*field[ind1+0]
						  + 3.0*field[ind1-1])*p1/sigma_rm2;

      ind1 = n1r-1;
      dfield[ind1] =  (field[ind1+3]
						  - 6.0*field[ind1+2]
						  + 15.0*field[ind1+1]
						  - 19.0*field[ind1+0]
						  + 12.0*field[ind1-1]
						  - 3.0*field[ind1-2])*p1/sigma_rm1;

      // right one
      ind1 = n1s+1;
      dfield[ind1] = (field[ind1-3]
						  - 6.0*field[ind1-2]
						  + 15.0*field[ind1-1]
						  - 19.0*field[ind1+0]
						  + 12.0*field[ind1+1]
						  - 3.0*field[ind1+2])*p1/sigma_sp1;

      ind1 = n1s+2;
      dfield[ind1] = (field[ind1-3]
						  - 6.0*field[ind1-2]
						  + 12.0*field[ind1-1]
						  - 10.0*field[ind1+0]
						  + 3.0*field[ind1+1])*p1/sigma_sp2;


      ind1 = n1s+3;
      dfield[ind1] =  (field[ind1-3]
						  - 3.0*field[ind1-2]
						  + 3.0*field[ind1-1]
						  - field[ind1+0])*p1/sigma_sp3;

    }

}





//-------------------------- KO-type dissipation, eight derivative in the interior --------------------------//
// Written by Manuel.

void diss_KO8_1d(struct GRID_PAR *grid,
	      FLOAT *field,
		 FLOAT *dfield)
{// diss_KO8

#ifdef EXCISION
  printf("************ THIS DIFFERENCE OPERATOR DOES NOT WORK WITH EXCISION !!! **************");
#endif



  int n1i = (*grid).start_grid;
  int n1f = (*grid).final_grid;

  FLOAT xi = (*grid).initial_x;
  FLOAT xf = (*grid).final_x;

  FLOAT p1;

  FLOAT h_1 =(xf-xi)/ (FLOAT)(n1f-n1i-1); // this is dx

  int n1r;
  int n1s;

  // scalar product
  FLOAT sigma_rm4;
  FLOAT sigma_rm3;
  FLOAT sigma_rm2;
  FLOAT sigma_rm1;
  FLOAT sigma_sp4;
  FLOAT sigma_sp3;
  FLOAT sigma_sp2;
  FLOAT sigma_sp1;

  // define what scalar product, interior means (ie which are gridpoints r,s), depending on the derivatives that we are using

  char macro_value_strg[100];
  GET_MACRO_VALUE(DERIV);


  if (strcmp(macro_value_strg,"derivS_1d")==0) { // second order in the interior

    n1r = n1i+4;  n1s = n1f-5;
    p1 = 1.0/pow(h_1,5);

    sigma_rm4 = 0.5; sigma_rm3 = 1.0; sigma_rm2 = 1.0;  sigma_rm1 = 1.0;
    sigma_sp4 = sigma_rm4;  sigma_sp3 = sigma_rm3;
    sigma_sp2 = sigma_rm2;  sigma_sp1 = sigma_rm1;
  }

  else if (strcmp(macro_value_strg,"derivQ_1d")==0) {// fourth order in the interior

    n1r = n1i+4;  n1s = n1f-5;
    p1 = 1.0/pow(h_1,3);

    sigma_rm4 = 17.0/48.0; sigma_rm3 = 59.0/48.0; sigma_rm2 = 43.0/48.0;  sigma_rm1 = 49.0/48.0;
    sigma_sp4 = sigma_rm4;  sigma_sp3 = sigma_rm3;
    sigma_sp2 = sigma_rm2;  sigma_sp1 = sigma_rm1;
  }

  else if (strcmp(macro_value_strg,"deriv_strand_third_order_boundaries_sixth_interior_1d")==0) { // sixth order in the interior

    n1r = n1i+6;  n1s = n1f-7;
    p1 = 1.0/h_1;

    sigma_rm4 = 2711.0/4320.0; sigma_rm3 = 5359.0/4320.0;
    sigma_rm2 = 7877.0/8640.0;  sigma_rm1 = 43801.0/43200.0;

    sigma_sp4 = sigma_rm4;  sigma_sp3 = sigma_rm3;
    sigma_sp2 = sigma_rm2;  sigma_sp1 = sigma_rm1;
  }

  else if (strcmp(macro_value_strg,"deriv_strand_fourth_order_boundaries_eight_interior_1d")==0) {// eighth order in the interior

   n1r = n1i+8;  n1s = n1f-9;
   p1 = h_1;

  sigma_rm4 = 299527.0/725760.0; sigma_rm3 = 103097.0/80640.0;
  sigma_rm2 = 670091.0/725760.0;  sigma_rm1 = 5127739.0/5080320.0;

  sigma_sp4 = sigma_rm4;  sigma_sp3 = sigma_rm3;
  sigma_sp2 = sigma_rm2;  sigma_sp1 = sigma_rm1;
  }

  else {
    printf("I don't know what derivative you are using, will stop here.\n");
    exit(1);
  }


  // initialize to zero
  {register int ind1;
  #pragma omp parallel for
  for (ind1 = n1i; ind1 < n1f-1; ind1++){

	dfield[ind1] =  0.0;

      } }


  //------------------------- dissipation in the x-direction ----------------------------//
  // interior
  {register int ind1;
  #pragma omp parallel for
  for (ind1 = n1r; ind1 < n1s+1; ind1++){
	dfield[ind1] = -(        field[ind1+4]
				    -  8.0*field[ind1+3]
				    + 28.0*field[ind1+2]
				    - 56.0*field[ind1+1]
				    + 70.0*field[ind1+0]
				    - 56.0*field[ind1-1]
				    + 28.0*field[ind1-2]
				    -  8.0*field[ind1-3]
				    +      field[ind1-4]
				    )*p1;
      } }

  // boundaries
  {register int ind1;

      // left one
      ind1 = n1r-4;
      dfield[ind1] =  ( -     field[ind1+4]
				  + 4.0*field[ind1+3]
				  - 6.0*field[ind1+2]
				  + 4.0*field[ind1+1]
				  -     field[ind1+0])*p1/sigma_rm4;
      ind1 = n1r-3;
      dfield[ind1] =   (    2.0*field[ind1+3]
				   -  9.0*field[ind1+2]
				   + 15.0*field[ind1+1]
				   - 11.0*field[ind1+0]
				   +  3.0*field[ind1-1])*p1/sigma_rm3;

      ind1 = n1r-2;
      dfield[ind1] =   ( -     field[ind1+3]
				   + 3.0*field[ind1+2]
				   - 8.0*field[ind1+0]
				   + 9.0*field[ind1-1]
				   - 3.0*field[ind1-2])*p1/sigma_rm2;

      ind1 = n1r-1;
      dfield[ind1] =  ( -  1.0*field[ind1+3]
				  +  6.0*field[ind1+2]
				  - 14.0*field[ind1+1]
				  + 15.0*field[ind1+0]
				  -  6.0*field[ind1-1]
				  -  1.0*field[ind1-2]
				  +  1.0*field[ind1-3])*p1/sigma_rm1;

      // right one
      ind1 = n1s+4;
      dfield[ind1] =  ( -field[ind1-4]
				  + 4.0*field[ind1-3]
				  - 6.0*field[ind1-2]
				  + 4.0*field[ind1-1]
				  - field[ind1+0])*p1/sigma_sp4;
      ind1 = n1s+3;
      dfield[ind1] =   ( 2.0*field[ind1-3]
				   - 9.0*field[ind1-2]
				   + 15.0*field[ind1-1]
				   - 11.0*field[ind1+0]
				   + 3.0*field[ind1+1])*p1/sigma_sp3;

      ind1 = n1s+2;
      dfield[ind1] =   ( -field[ind1-3]
				   + 3.0*field[ind1-2]
				   - 8.0*field[ind1+0]
				   + 9.0*field[ind1+1]
				   - 3.0*field[ind1+2])*p1/sigma_sp2;

      ind1 = n1s+1;
      dfield[ind1] =  ( - 1.0*field[ind1-3]
				  + 6.0*field[ind1-2]
				  - 14.0*field[ind1-1]
				  + 15.0*field[ind1+0]
				  - 6.0*field[ind1+1]
				  - 1.0*field[ind1+2]
				  + 1.0*field[ind1+3])*p1/sigma_sp1;

    }

 }





/* ------------------------------------------------------------------------------------------ */



/* ------------------------------------------------------------------------------------------ */

//-------------------------- KO-type dissipation, fourth derivative Periodic --------------------------//
// Written by Manuel.
void diss_KO4_Per_1d(struct GRID_PAR *grid,
	      FLOAT *field,
		  FLOAT *dfield)
{// diss_KO4

#ifdef EXCISION
  printf("************ THIS DIFFERENCE OPERATOR DOES NOT WORK WITH EXCISION !!! **************");
#endif


  int n1i = (*grid).start_grid;
  int n1f = (*grid).final_grid;
  int n_tot = n1f - n1i;

  FLOAT xi = (*grid).initial_x;
  FLOAT xf = (*grid).final_x;

  FLOAT d2 = 2.0;
  FLOAT d4 = 4.0;
  FLOAT d5 = 5.0;
  FLOAT d6 = 6.0;

  FLOAT p1;


  FLOAT h_1 = (xf-xi)/(FLOAT)(n1f-n1i); // this is dx


  // define what scalar product, interior means (ie which are gridpoints r,s), depending on the derivatives that we are using

  char macro_value_strg[100];
  GET_MACRO_VALUE(DERIV);

  // --------- derivatives on the sphere ------------
  if (strcmp(macro_value_strg,"derivD_Per_1d")==0) { // second order in the interior
    p1 = 1.0/h_1;
  }
  else if (strcmp(macro_value_strg,"derivQ_Per_1d")==0) {// fourth order in the interior
    p1 = h_1;
  }
  else if (strcmp(macro_value_strg,"deriv_6_Per_1d")==0) { // sixth order in the interior
    p1 = pow(h_1,3);
  }
  else if (strcmp(macro_value_strg,"deriv_8_Per_1d")==0) {// eighth order in the interior
  p1 = pow(h_1,5);
  }
  else {
    printf("I don't know what derivative you are using, will stop here.\n");
    exit(1);
  }


  // initialize to zero

  {register int ind1;
  for (ind1 = n1i; ind1 < n1f; ind1++){
	dfield[ind1] =  0.0;
      }
    }


  //------------------------- dissipation in the x-direction ----------------------------//
  // interior
  {register int ind1;
  	  #ifdef OMP_D
      #pragma omp parallel for  default(shared) private(ind1)
      #endif //OMP_D
  for (ind1 = n1i; ind1 < n1f; ind1++){

	dfield[ind1] =  -(     field[(ind1-2 + n_tot) % n_tot]
								- d4*field[(ind1-1 + n_tot) % n_tot]
								+ d6*field[ind1+0]
								- d4*field[(ind1+1) % n_tot]
								+    field[(ind1+2) % n_tot]
								)*p1;
	//printf("i-2=%d,i-1=%d,i=%d,i+1=%d,i+2=%d \n",)
      }
      #ifdef OMP_D

      #endif //OMP_D
      }





}


//-------------------------- KO-type dissipation, sixth derivative in the interior --------------------------//
// Written by Manuel.

void diss_KO6_Per_1d(struct GRID_PAR *grid,
	      FLOAT *field,
		  FLOAT *dfield)
{// diss_KO6




  int n1i = (*grid).start_grid;
  int n1f = (*grid).final_grid;
  int n_tot = n1f - n1i;

  FLOAT xi = (*grid).initial_x;
  FLOAT xf = (*grid).final_x;




  FLOAT p1;

  FLOAT h_1 = (xf-xi)/(FLOAT)(n1f-n1i); // this is dx





  // define what scalar product, interior means (ie which are gridpoints r,s), depending on the derivatives that we are using

  char macro_value_strg[100];
  GET_MACRO_VALUE(DERIV);

  // --------- derivatives on the sphere ------------
  if (strcmp(macro_value_strg,"derivD_Per_1d")==0) { // second order in the interior
    p1 = 1.0/pow(h_1,3);
  }

  else if (strcmp(macro_value_strg,"derivQ_Per_1d")==0) {// fourth order in the interior
    p1 = 1.0/h_1;
  }

  else if (strcmp(macro_value_strg,"deriv_6_Per_1d")==0) { // sixth order in the interior
    p1 = h_1;
  }

  else if (strcmp(macro_value_strg,"deriv_8_Per_1d")==0) {// eighth order in the interior
  p1 = pow(h_1,3);
  }

  else {
    printf("I don't know what derivative you are using, will stop here.\n");
    exit(1);
  }

#ifdef EXCISION
  printf("************ THIS DIFFERENCE OPERATOR DOES NOT WORK WITH EXCISION !!! **************");
#endif

  // initialize to zero
  {register int ind1;
  for (ind1 = n1i; ind1 < n1f-1; ind1++){

	dfield[ind1] =  0.0;

      } }


  //------------------------- dissipation in the x-direction ----------------------------//
  // interior
  {register int ind1;
      #ifdef OMP_D
      #pragma omp parallel for  default(shared) private(ind1)
      #endif //OMP_D
  for (ind1 = n1i; ind1 < n1f; ind1++){
	dfield[ind1] =  ( field[(ind1+3 + n_tot) % n_tot]
						     - 6.0*field[(ind1+2 + n_tot) % n_tot]
						     + 15.0*field[(ind1+1 + n_tot) % n_tot]
						     - 20.0*field[ind1+0]
						     + 15.0*field[(ind1-1 + n_tot) % n_tot]
						     - 6.0*field[(ind1-2 + n_tot) % n_tot]
						     + field[(ind1-3 + n_tot) % n_tot])*p1;
      }
      #ifdef OMP_D

      #endif //OMP_D
      }

}





//-------------------------- KO-type dissipation, eight derivative in the interior --------------------------//
// Written by Manuel.

void diss_KO8_Per_1d(struct GRID_PAR *grid,
	      FLOAT *field,
		  FLOAT *dfield)
{// diss_KO8

#ifdef EXCISION
  printf("************ THIS DIFFERENCE OPERATOR DOES NOT WORK WITH EXCISION !!! **************");
#endif



  int n1i = (*grid).start_grid;
  int n1f = (*grid).final_grid;

  FLOAT xi = (*grid).initial_x;
  FLOAT xf = (*grid).final_x;
  int n_tot = n1f - n1i;

  FLOAT p1;

  FLOAT h_1 =(xf-xi)/ (FLOAT)(n1f-n1i); // this is dx

  // define what scalar product, interior means (ie which are gridpoints r,s), depending on the derivatives that we are using

  char macro_value_strg[100];
  GET_MACRO_VALUE(DERIV);

  // --------- derivatives on the sphere ------------
  if (strcmp(macro_value_strg,"derivD_Per_1d")==0) { // second order in the interior
    p1 = 1.0/pow(h_1,5);
  }

  else if (strcmp(macro_value_strg,"derivQ_Per_1d")==0) {// fourth order in the interior
    p1 = 1.0/pow(h_1,3);
  }

  else if (strcmp(macro_value_strg,"deriv_6_Per_1d")==0) { // sixth order in the interior
    p1 = 1.0/h_1;
  }

  else if (strcmp(macro_value_strg,"deriv_8_Per_1d")==0) {// eighth order in the interior
   p1 = h_1;
  }

  else {
    printf("I don't know what derivative you are using, will stop here.\n");
    exit(1);
  }


  // initialize to zero
  {register int ind1;
      #ifdef OMP_D
      #pragma omp parallel for  default(shared) private(ind1)
      #endif //OMP_D
  for (ind1 = n1i; ind1 < n1f-1; ind1++){

	dfield[ind1] =  0.0;

      }
      #ifdef OMP_D

      #endif //OMP_D
      }


  //------------------------- dissipation in the x-direction ----------------------------//
  // interior
  {register int ind1;
      #ifdef OMP_D
      #pragma omp parallel for  default(shared) private(ind1)
      #endif //OMP_D
  for (ind1 = n1i; ind1 < n1f; ind1++){
	dfield[ind1] =  -( field[(ind1+4 + n_tot) % n_tot]
						     - 8.0*field[(ind1+3 + n_tot) % n_tot]
						     + 28.0*field[(ind1+2 + n_tot) % n_tot]
						     - 56.0*field[(ind1+1 + n_tot) % n_tot]
						     + 70.0*field[ind1+0]
						     - 56.0*field[(ind1-1 + n_tot) % n_tot]
						     + 28.0*field[(ind1-2 + n_tot) % n_tot]
						     - 8.0*field[(ind1-3 + n_tot) % n_tot]
						     + field[(ind1-4 + n_tot) % n_tot]
						     )*p1;
      }
      #ifdef OMP_D

      #endif //OMP_D
      }
 }





#endif
