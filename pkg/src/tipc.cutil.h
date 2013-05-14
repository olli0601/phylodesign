/** \file tipc.cutil.h
    \brief Main file that provides all functions that are callable from R.
*/

#ifndef TIPC_CUTIL_H_
#define TIPC_CUTIL_H_

#include <R.h>
#include <Rinternals.h>

#define NEW(x) (x *)malloc(sizeof(x))									/**< macro to allocate new memory*/
#define NEW_ARY(x,y) (x *)malloc((y)*sizeof(x))					/**< macro to allocate a contiguous memory array*/
#define NEW_ZERO_ARY(x,y) (x *)calloc(y,sizeof(x))			/**< macro to allocate a contiguous memory array that is initialized to zero*/
#define RLC_ARY(x,y,z) x = (z *)realloc((x),(y)*sizeof(z))		/**< macro to re-allocate a contiguous memory array; this preserves the content while resizing*/
#define DELETE(x) free(x)														/**< macro to free memory*/
#define MOVE(ty, fr, to, n) (ty *)memmove(to, fr, (n)*sizeof(ty))	/** < macro to move an array */	//void * memmove (void *to, const void *from, size_t size)
#define COPY_INTO_SEPARATE(ty, fr, to, n)	(ty *)memcpy(to, fr, (n)*sizeof(ty))		/** < macro to copy an array into a non-overlapping array; both arrays are assumed to be at least 'n' units long */			//void * memcpy ( void * destination, const void * source, size_t num );
#define CAST(x,y) (static_cast<x>(y))									/**< macro to cast between types*/
#define MIN(x,y) ((x < y) ? x : y) 											/**< macro to return the minimum of two values */
#define MAX(x,y) ((x < y) ? y : x)											/**< macro to return the maximum of two values */
#define ABSDIFF(x,y) ((y-x>0) ? y-x : x-y) 								/**< macro to return the absolute difference of two values */
#define ABS(x) ((x>0) ? x : -x) 												/**< macro to return the absolute value of a number */


extern "C" {

SEXP tipc_tabulate_after_sample(SEXP tpc, SEXP sid);

}

#endif /*TIPC_CUTIL_H_*/
