/* binio.h */


/*
 * Functions to do binary I/O of floats, ints, etc. with byte swapping
 * as needed.
 */

#ifndef BINIO_H
#define BINIO_H

/* Include files which define SEEK_SET, O_RD_ONLY, etc. */
/* and prototype open(), close(), lseek(), etc. */
#include <unistd.h>
#include <fcntl.h>

#ifdef __STDC__
#  define PROTO(x) x
#else
#  define PROTO(x) ()	
#endif

extern void flip4 PROTO(( unsigned int *src, unsigned int *dest, int n ));
extern void flip2 PROTO(( unsigned short *src, unsigned short *dest, int n ));

#ifdef cray
  extern void cray_to_ieee_array( /* dest, source, n */ );
  extern void ieee_to_cray_array( /* dest, source, n */ );
#endif

/******************************************************/
/*                    Read Functions                  */
/******************************************************/
extern int read_int4 PROTO(( int f, int *i ));
extern int read_float4 PROTO(( int f, float *x ));

/*******************************************************/
/*               Write Functions                       */
/*******************************************************/
extern int write_int4 PROTO(( int f, int i ));
extern int write_float4 PROTO(( int f, float x ));

#undef PROTO

#endif
