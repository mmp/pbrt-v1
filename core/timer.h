
/*
 * pbrt source code Copyright(c) 1998-2005 Matt Pharr and Greg Humphreys
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 * (See file License.txt for complete license)
 */

#ifndef PBRT_TIMER_H
#define PBRT_TIMER_H
// timer.h*
#include "pbrt.h"
#if defined ( WIN32 )
#include <windows.h>
#elif defined (IRIX) || defined (IRIX64)
#include <stddef.h>
#include <fcntl.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <sys/syssgi.h>
#include <sys/errno.h>
#include <unistd.h>
#else
#include <sys/time.h>
#endif
// Timer Declarations
class COREDLL Timer {
public:
	// Public Timer Methods
	Timer();
	~Timer();
	
	void Start();
	void Stop();
	void Reset();
	
	double Time();
private:
	// Private Timer Data
	double time0, elapsed;
	bool running;
	double GetTime();
#if defined( IRIX ) || defined( IRIX64 )
	// Private IRIX Timer Data
	int fd;
	unsigned long long counter64;
	unsigned int counter32;
	unsigned int cycleval;
	
	typedef unsigned long long iotimer64_t;
	typedef unsigned int iotimer32_t;
	volatile iotimer64_t *iotimer_addr64;
	volatile iotimer32_t *iotimer_addr32;
	
	void *unmapLocation;
	int unmapSize;
#elif defined( WIN32 )
	// Private Windows Timer Data
	LARGE_INTEGER performance_counter, performance_frequency;
	double one_over_frequency;
#else
	// Private UNIX Timer Data
	struct timeval timeofday;
#endif
};
#endif // PBRT_TIMER_H
