#if defined(SGI)
/* 
 --- Fortran-callable routine ZUNDER that sets the bit to specify
 --- that underflows are flushed to zero in hardware on SGI R10000.
 --- See man handle_sigfpes
 --- Alan J. Wallcraft,  NRL,  October 1997.
*/
#if defined(NERSC_VERSION)
/*Code not suported on Ve and Vilje*/
#else
#include <sys/fpu.h>
void zunder_()
{
	union fpc_csr   n;
	n.fc_word = get_fpc_csr();
	n.fc_struct.flush = 1;
	set_fpc_csr(n.fc_word);
}
#endif /* NERSC_VERSION */
#endif /* SGI */

#if defined(AIX)
/* 
 --- Fortran-callable function WTIME that returns the wall time in seconds.
 --- Probably not thread-safe, only for Power-PC systems.
 --- Alan J. Wallcraft,  NRL,  May 2001.
 --- Based on notes by Bob Walkup (10x faster than MPI_WTIME).
*/
#include <sys/time.h>
#include <sys/systemcfg.h>
double wtime(void)
{
	struct timebasestruct TB;
	static int    first_call;
	static double tb_factor;
	double tb_top,tb_bot;
	if (first_call == 0) {
		first_call = 1;
		tb_top     = (double) _system_configuration.Xint;
		tb_bot     = (double) _system_configuration.Xfrac;
		tb_factor  = tb_top/tb_bot;
        }
	read_real_time(&TB, TIMEBASE_SZ);
	return ( tb_factor * ( 4.294967296*((double) TB.tb_high) + 1.0e-9*((double) TB.tb_low) ) );
}
#endif /* AIX */

#if defined(ENDIAN_IO)
/* 
 --- Fortran-callable routine ZAIO_ENDIAN to swap the endian-ness of an array
 --- Brent Anderson, ASC MSRC, July 2007.
*/
#include <inttypes.h>
void zaio_endian_(uint32_t *, uint32_t *);

#define swap4bytes(data) \
 ( (((data) >> 24) & 0x000000FF) | \
   (((data) >>  8) & 0x0000FF00) | \
   (((data) <<  8) & 0x00FF0000) | \
   (((data) << 24) & 0xFF000000) )

/** Swaps the byte-order in a 32-bit word from
    big- to little-endianness or vice-versa. */

void zaio_endian_(uint32_t aa[], uint32_t *nn) {
    uint32_t ii;
    for( ii=0; ii<*nn; ii++ ) {
        aa[ii] = swap4bytes ( aa[ii] );
    }
    return;
}
#endif /* ENDIAN_IO */

void machine_c()
{
}
#if defined(NERSC_VERSION)
/*
 --- These are wrapper routines for libc system calls. Not necessary on some arch/compiler
 --- combos, but using these wrappers makes the code more portable
 --- NB: Note the underscore - It is added because fortran usually  appends an underscore 
 ----    after symbol names. If the linker complains about not finding symbols, try to 
 ----    remove the underscore
*/
#include <sys/stat.h>
#include <sys/types.h>
#if defined(AIX)
int mkdir_wrap (char *dirname)
#else
int mkdir_wrap_(char *dirname)
#endif
{
   if (mkdir(dirname, 0777) == 0) {
      return -1;
   } else {
      return 0;
   }
}


#include <unistd.h>
#if defined(AIX)
int unlink_wrap (char *fname)
#else
int unlink_wrap_(char *fname)
#endif
{
   if (unlink(fname) == 0) {
      return -1;
   } else {
      return 0;
   }
}
#endif /* NERSC_VERSION */
