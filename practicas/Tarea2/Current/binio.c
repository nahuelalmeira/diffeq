/* binio.c */

/*
 * The file format is assumed to be BIG-ENDIAN.
 * If this code is compiled with -DLITTLE and executes on a little endian
 * CPU then byte-swapping will be done.
 *
 * If an ANSI compiler is used prototypes and ANSI function declarations
 * are used.  Otherwise use K&R conventions.
 *
 * If we're running on a CRAY (8-byte ints and floats), conversions will
 * be done as needed.
 */

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#ifdef cray
#  include <string.h>
#endif
#include "binio.h"


/**********************************************************/
/*                  Byte Flipping                         */
/**********************************************************/


#define FLIP4( n )  (  (n & 0xff000000) >> 24		\
		     | (n & 0x00ff0000) >> 8		\
		     | (n & 0x0000ff00) << 8		\
		     | (n & 0x000000ff) << 24  )


/*
 * Flip the order of the 4 bytes in an array of 4-byte words.
 */
#ifdef __STDC__
  void flip4( unsigned int *src, unsigned int *dest, int n )
#else
  void flip4( src, dest, n )
  unsigned int *src, *dest;
  int n;
#endif
{
   int i;
   for (i=0;i<n;i++) {
      register unsigned int tmp = src[i];
      dest[i] = FLIP4( tmp );
   }
}


#ifdef cray

/** THESE ROUTINES MUST BE COMPILED ON THE CRAY ONLY SINCE THEY **/
/** REQUIRE 8-BYTES PER C-TYPE LONG                             **/

/* Cray to IEEE single precision */
static void c_to_if(t,f)
long *t,*f;
{
    if (*f != 0){
        *t = (((*f & 0x8000000000000000) |      /* sign bit */
                 ((((*f & 0x7fff000000000000) >> 48)-16258) << 55)) + /* exp */
                 (((*f & 0x00007fffff000000) +
                    ((*f & 0x0000000000800000) << 1)) << 8));  /* mantissa */
    }
    else *t = *f;
}


#define C_TO_IF( T, F )							\
	if (F != 0) {							\
		T = (((F & 0x8000000000000000) |			\
		((((F & 0x7fff000000000000) >> 48)-16258) << 55)) +	\
		(((F & 0x00007fffff000000) +				\
		((F & 0x0000000000800000) << 1)) << 8));		\
	}								\
	else {								\
		T = F;							\
	}



/* IEEE single precison to Cray */
static void if_to_c(t,f)
long *t,*f;
{
    if (*f != 0) {
        *t = (((*f & 0x8000000000000000) |
                ((*f & 0x7f80000000000000) >> 7) +
                (16258 << 48)) |
                (((*f & 0x007fffff00000000) >> 8) | (0x0000800000000000)));
        if ((*f << 1) == 0) *t = 0;
    }
    else *t = *f;
}

/* T and F must be longs! */
#define IF_TO_C( T, F )							\
	if (F != 0) {							\
		T = (((F & 0x8000000000000000) |			\
		((F & 0x7f80000000000000) >> 7) +			\
		(16258 << 48)) |					\
		(((F & 0x007fffff00000000) >> 8) | (0x0000800000000000)));  \
		if ((F << 1) == 0) T = 0;				\
	}								\
	else {								\
		T = F;							\
	}

#endif /*cray*/



/**************************************************************/
/*                     Read Functions                         */
/**************************************************************/

/*
 * receive_data 
 *  Receive a block of data by repeatedly reading from a socket.  This
 *  function will block the process if data is not available.
 *  Input:  socket - socket to read from
 *           buffer - pointer to data buffer
 *          bytes - how many bytes to read into buffer
 *  Return:  number of bytes received.  This will always either be
 *           the same as the bytes requested or 0 if the connection
 *           was lost.
 */
int receive_data( f, buffer, bytes )
int f;
char *buffer;
int bytes;
{
   int sofar, remaining, len;

   sofar = 0;
   remaining = bytes;
   do {
      len = read( f, buffer + sofar, remaining);
      if (len<=0) return 0;
      sofar += len;
      remaining -= len;
   } while (remaining > 0);
   return bytes;
}


/*
 * Read a 4-byte integer. This function will
 *  block if the data is not available.
 * Input:  f - the file descriptor to read from
 *         i - pointer to integer to put result into.
 * Return:  1 = ok, 0 = error
 */
#ifdef __STDC__
  int read_int4( int f, int *i )
#else
  int read_int4( f, i, 4 )
  int f;
  int *i;
#endif
{
#ifdef LITTLE
   /* read big endian and convert to little endian */
   unsigned int n;
   if (receive_data( f, &n, 4 )==4) {
      *i = FLIP4( n );
      return 1;
   }
   else {
      return 0;
   }
#else
   if (receive_data( f, i, 4 )==4) {
#  ifdef cray
      *i = *i >> 32;
#  endif
      return 1;
   }
   else {
      return 0;
   }
#endif
}


/*
 * Read a 4-byte IEEE float.
 * Input:  f - the file descriptor to read from.
 *         x - pointer to float to put result into.
 * Return:  1 = ok, 0 = error
 */
#ifdef __STDC__
  int read_float4( int f, float *x )
#else
  int read_float4( f, x )
  int f;
  float *x;
#endif
{
#ifdef cray
   long buffer = 0;

   if ( receive_data( f, &buffer, 4 )==4 ) {
      /* convert IEEE float (buffer) to Cray float (x) */
      if_to_c( x, &buffer );
      return 1;
    }
    return 0;
#else
#  ifdef LITTLE
      unsigned int n, *iptr;
      if (receive_data( f, &n, 4 )==4) {
	 iptr = (unsigned int *) x;
	 *iptr = FLIP4( n );
	 return 1;
      }
      else {
	 return 0;
      }
#  else
      if (receive_data( f, x, 4 )==4) {
	 return 1;
      }
      else {
	 return 0;
      }
#  endif
#endif
}


/********************************************************/
/*                         Write Functions              */
/********************************************************/

/*
 * Write a 4-byte integer.
 *Input:  f - the file descriptor
 *         i - the integer
 * Return:  1 = ok, 0 = error
 */
#ifdef __STDC__
  int write_int4( int f, int i )
#else
  int write_int4( f, i )
  int f;
  int i;
#endif
{
#ifdef cray
   i = i << 32;
   return write( f, &i, 4 ) > 0;
#else
#  ifdef LITTLE
     i = FLIP4( i );
#  endif
   return write( f, &i, 4 ) > 0;
#endif
}


/*
 * Write a 4-byte IEEE float.
 * Input:  f - the file descriptor
 *         x - the float
 * Return:  1 = ok, 0 = error
 */
#ifdef __STDC__
  int write_float4( int f, float x )
#else
  int write_float4( f, x )
  int f;
  float x;
#endif
{
#ifdef cray
   char buffer[8];
   c_to_if( buffer, &x );
   return write( f, buffer, 4 ) > 0;
#else
#  ifdef LITTLE
      float y;
      unsigned int *iptr = (unsigned int *) &y, temp;
      y = (float) x;
      temp = FLIP4( *iptr );
      return write( f, &temp, 4 ) > 0;
#  else
      float y;
      y = (float) x;
      return write( f, &y, 4 ) > 0;
#  endif
#endif
}

