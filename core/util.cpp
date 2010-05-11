
/*
    pbrt source code Copyright(c) 1998-2010 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    pbrt is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.  Note that the text contents of
    the book "Physically Based Rendering" are *not* licensed under the
    GNU GPL.

    pbrt is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

 */

// util.cpp*
#include "pbrt.h"
#include "timer.h"
#include <map>
using std::map;
// Error Reporting Includes
#include <stdarg.h>
// Error Reporting Definitions
#define PBRT_ERROR_IGNORE 0
#define PBRT_ERROR_CONTINUE 1
#define PBRT_ERROR_ABORT 2
// Error Reporting Functions
static void processError(const char *format, va_list args,
		const char *message, int disposition) {
#ifndef WIN32
	char *errorBuf;
	vasprintf(&errorBuf, format, args);
#else
	char errorBuf[2048];
	_vsnprintf(errorBuf, sizeof(errorBuf), format, args);
#endif
	// Report error
	switch (disposition) {
	case PBRT_ERROR_IGNORE:
		return;
	case PBRT_ERROR_CONTINUE:
		fprintf(stderr, "%s: %s\n", message, errorBuf);
		// Print scene file and line number, if appropriate
		extern int line_num;
		if (line_num != 0) {
			extern string current_file;
			fprintf(stderr, "\tLine %d, file %s\n", line_num,
				current_file.c_str());
		}
		break;
	case PBRT_ERROR_ABORT:
		fprintf(stderr, "%s: %s\n", message, errorBuf);
		// Print scene file and line number, if appropriate
		extern int line_num;
		if (line_num != 0) {
			extern string current_file;
			fprintf(stderr, "\tLine %d, file %s\n", line_num,
				current_file.c_str());
		}
		abort();
	}
#ifndef WIN32
	free(errorBuf);
#endif
}
COREDLL void Info(const char *format, ...) {
	va_list args;
	va_start(args, format);
	processError(format, args, "Notice", PBRT_ERROR_CONTINUE);
	va_end(args);
}
COREDLL void Warning(const char *format, ...) {
	va_list args;
	va_start(args, format);
	processError(format, args, "Warning", PBRT_ERROR_CONTINUE);
	va_end(args);
}
COREDLL void Error(const char *format, ...) {
	va_list args;
	va_start(args, format);
	processError(format, args, "Error", PBRT_ERROR_CONTINUE);
	va_end(args);
}
COREDLL void Severe(const char *format, ...) {
	va_list args;
	va_start(args, format);
	processError(format, args, "Fatal Error", PBRT_ERROR_ABORT);
	va_end(args);
}
// Matrix Method Definitions
COREDLL bool SolveLinearSystem2x2(const float A[2][2],
		const float B[2], float x[2]) {
	float det = A[0][0]*A[1][1] - A[0][1]*A[1][0];
	if (fabsf(det) < 1e-5)
		return false;
	float invDet = 1.0f/det;
	x[0] = (A[1][1]*B[0] - A[0][1]*B[1]) * invDet;
	x[1] = (A[0][0]*B[1] - A[1][0]*B[0]) * invDet;
	return true;
}
Matrix4x4::Matrix4x4(float mat[4][4]) {
	memcpy(m, mat, 16*sizeof(float));
}
Matrix4x4::Matrix4x4(float t00, float t01, float t02, float t03,
                     float t10, float t11, float t12, float t13,
                     float t20, float t21, float t22, float t23,
                     float t30, float t31, float t32, float t33) {
	m[0][0] = t00; m[0][1] = t01; m[0][2] = t02; m[0][3] = t03;
	m[1][0] = t10; m[1][1] = t11; m[1][2] = t12; m[1][3] = t13;
	m[2][0] = t20; m[2][1] = t21; m[2][2] = t22; m[2][3] = t23;
	m[3][0] = t30; m[3][1] = t31; m[3][2] = t32; m[3][3] = t33;
}
Reference<Matrix4x4> Matrix4x4::Transpose() const {
   return new Matrix4x4(m[0][0], m[1][0], m[2][0], m[3][0],
	                    m[0][1], m[1][1], m[2][1], m[3][1],
	                    m[0][2], m[1][2], m[2][2], m[3][2],
	                    m[0][3], m[1][3], m[2][3], m[3][3]);
}
Reference<Matrix4x4> Matrix4x4::Inverse() const {
	int indxc[4], indxr[4];
	int ipiv[4] = { 0, 0, 0, 0 };
	float minv[4][4];
	memcpy(minv, m, 4*4*sizeof(float));
	for (int i = 0; i < 4; i++) {
		int irow = -1, icol = -1;
		float big = 0.;
		// Choose pivot
		for (int j = 0; j < 4; j++) {
			if (ipiv[j] != 1) {
				for (int k = 0; k < 4; k++) {
					if (ipiv[k] == 0) {
						if (fabsf(minv[j][k]) >= big) {
							big = float(fabsf(minv[j][k]));
							irow = j;
							icol = k;
						}
					}
					else if (ipiv[k] > 1)
						Error("Singular matrix in MatrixInvert");
				}
			}
		}
		++ipiv[icol];
		// Swap rows _irow_ and _icol_ for pivot
		if (irow != icol) {
			for (int k = 0; k < 4; ++k)
				swap(minv[irow][k], minv[icol][k]);
		}
		indxr[i] = irow;
		indxc[i] = icol;
		if (minv[icol][icol] == 0.)
			Error("Singular matrix in MatrixInvert");
		// Set $m[icol][icol]$ to one by scaling row _icol_ appropriately
		float pivinv = 1.f / minv[icol][icol];
		minv[icol][icol] = 1.f;
		for (int j = 0; j < 4; j++)
			minv[icol][j] *= pivinv;
		// Subtract this row from others to zero out their columns
		for (int j = 0; j < 4; j++) {
			if (j != icol) {
				float save = minv[j][icol];
				minv[j][icol] = 0;
				for (int k = 0; k < 4; k++)
					minv[j][k] -= minv[icol][k]*save;
			}
		}
	}
	// Swap columns to reflect permutation
	for (int j = 3; j >= 0; j--) {
		if (indxr[j] != indxc[j]) {
			for (int k = 0; k < 4; k++)
				swap(minv[k][indxr[j]], minv[k][indxc[j]]);
		}
	}
	return new Matrix4x4(minv);
}
// Statistics Definitions
struct COREDLL StatTracker {
	StatTracker(const string &cat, const string &n,
	            StatsCounterType *pa, StatsCounterType *pb = NULL,
		    bool percentage = true);
	string category, name;
	StatsCounterType *ptra, *ptrb;
	bool percentage;
};
typedef map<std::pair<string, string>, StatTracker *> TrackerMap;
static TrackerMap trackers;
static void addTracker(StatTracker *newTracker) {
	std::pair<string, string> s = std::make_pair(newTracker->category, newTracker->name);
	if (trackers.find(s) != trackers.end()) {
		newTracker->ptra = trackers[s]->ptra;
		newTracker->ptrb = trackers[s]->ptrb;
		return;
	}
	trackers[s] = newTracker;
}
static void StatsPrintVal(FILE *f, StatsCounterType v);
static void StatsPrintVal(FILE *f, StatsCounterType v1, StatsCounterType v2);
// Statistics Functions
StatTracker::StatTracker(const string &cat, const string &n,
                         StatsCounterType *pa, StatsCounterType *pb, bool p) {
	category = cat;
	name = n;
	ptra = pa;
	ptrb = pb;
	percentage = p;
}
StatsCounter::StatsCounter(const string &category, const string &name) {
	num = 0;
	addTracker(new StatTracker(category, name, &num));
}
StatsRatio::StatsRatio(const string &category, const string &name) {
	na = nb = 0;
	addTracker(new StatTracker(category, name, &na, &nb, false));
}
StatsPercentage::StatsPercentage(const string &category, const string &name) {
	na = nb = 0;
	addTracker(new StatTracker(category, name, &na, &nb, true));
}
void StatsPrint(FILE *dest) {
	fprintf(dest, "Statistics:\n");
	TrackerMap::iterator iter = trackers.begin();
	string lastCategory;
	while (iter != trackers.end()) {
		// Print statistic
		StatTracker *tr = iter->second;
		if (tr->category != lastCategory) {
			fprintf(dest, "%s\n", tr->category.c_str());
			lastCategory = tr->category;
		}
		fprintf(dest, "    %s", tr->name.c_str());
		// Pad out to results column
		int resultsColumn = 56;
		int paddingSpaces = resultsColumn - (int) tr->name.size();
		while (paddingSpaces-- > 0)
			putc(' ', dest);
		if (tr->ptrb == NULL)
			StatsPrintVal(dest, *tr->ptra);
		else {
			if (*tr->ptrb > 0) {
				float ratio = (float)*tr->ptra / (float)*tr->ptrb;
				StatsPrintVal(dest, *tr->ptra, *tr->ptrb);
				if (tr->percentage)
					fprintf(dest, " (%3.2f%%)", 100. * ratio);
				else
					fprintf(dest, " (%.2fx)", ratio);
			}
			else
				StatsPrintVal(dest, *tr->ptra, *tr->ptrb);
		}
		fprintf(dest, "\n");
		++iter;
	}
}
static void StatsPrintVal(FILE *f, StatsCounterType v) {
	if (v > 1e9) fprintf(f, "%.3fB", v / 1e9f);
	else if (v > 1e6) fprintf(f, "%.3fM", v / 1e6f);
	else if (v > 1e4) fprintf(f, "%.1fk", v / 1e3f);
	else fprintf(f, "%.0f", (float)v);
}
static void StatsPrintVal(FILE *f, StatsCounterType v1,
		StatsCounterType v2) {
	StatsCounterType m = min(v1, v2);
	if (m > 1e9) fprintf(f, "%.3fB:%.3fB", v1 / 1e9f, v2 / 1e9f);
	else if (m > 1e6) fprintf(f, "%.3fM:%.3fM", v1 / 1e6f, v2 / 1e6f);
	else if (m > 1e4) fprintf(f, "%.1fk:%.1fk", v1 / 1e3f, v2 / 1e3f);
	else fprintf(f, "%.0f:%.0f", v1, v2);
}
void StatsCleanup() {
	TrackerMap::iterator iter = trackers.begin();
	string lastCategory;
	while (iter != trackers.end()) {
		delete iter->second;
		++iter;
	}
	trackers.erase(trackers.begin(), trackers.end());
}
// Random Number State
/*
   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote
        products derived from this software without specific prior written
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */
// Random Number Functions
static void init_genrand(u_long seed) {
	mt[0]= seed & 0xffffffffUL;
	for (mti=1; mti<N; mti++) {
		mt[mti] =
		(1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
		/* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
		/* In the previous versions, MSBs of the seed affect   */
		/* only MSBs of the array mt[].                        */
		/* 2002/01/09 modified by Makoto Matsumoto             */
		mt[mti] &= 0xffffffffUL;
		/* for >32 bit machines */
	}
}
COREDLL unsigned long genrand_int32(void)
{
	unsigned long y;
	static unsigned long mag01[2]={0x0UL, MATRIX_A};
	/* mag01[x] = x * MATRIX_A  for x=0,1 */

	if (mti >= N) { /* generate N words at one time */
		int kk;

		if (mti == N+1)   /* if init_genrand() has not been called, */
			init_genrand(5489UL); /* default initial seed */

		for (kk=0;kk<N-M;kk++) {
			y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
			mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
		}
		for (;kk<N-1;kk++) {
			y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
			mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
		}
		y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
		mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

		mti = 0;
	}

	y = mt[mti++];

	/* Tempering */
	y ^= (y >> 11);
	y ^= (y << 7) & 0x9d2c5680UL;
	y ^= (y << 15) & 0xefc60000UL;
	y ^= (y >> 18);

	return y;
}
/* generates a random number on [0,1)-real-interval */
COREDLL float genrand_real2(void)
{
	return (RandomUInt() & 0xffffff) / float(1 << 24);
}
// Memory Allocation Functions
COREDLL void *AllocAligned(size_t size) {
#ifndef L1_CACHE_LINE_SIZE
#define L1_CACHE_LINE_SIZE 64
#endif
	return memalign(L1_CACHE_LINE_SIZE, size);
}
COREDLL void FreeAligned(void *ptr) {
#ifdef WIN32 // NOBOOK
	_aligned_free(ptr);
#else // NOBOOK
	free(ptr);
#endif // NOBOOK
}
// ProgressReporter Method Definitions
ProgressReporter::ProgressReporter(int totalWork, const string &title, int bar_length)
	: totalPlusses(bar_length - title.size()) {
	plussesPrinted = 0;
	frequency = (float)totalWork / (float)totalPlusses;
	count = frequency;
	timer = new Timer();
	timer->Start();
	outFile = stdout;
	// Initialize progress string
	const int bufLen = title.size() + totalPlusses + 64;
	buf = new char[bufLen];
	snprintf(buf, bufLen, "\r%s: [", title.c_str());
	curSpace = buf + strlen(buf);
	char *s = curSpace;
	for (int i = 0; i < totalPlusses; ++i)
		*s++ = ' ';
	*s++ = ']';
	*s++ = ' ';
	*s++ = '\0';
	fprintf(outFile, buf);
	fflush(outFile);
}
ProgressReporter::~ProgressReporter() { delete[] buf; delete timer; }
void ProgressReporter::Update(int num) const {
	count -= num;
	bool updatedAny = false;
	while (count <= 0) {
		count += frequency;
		if (plussesPrinted++ < totalPlusses)
			*curSpace++ = '+';
		updatedAny = true;
	}
	if (updatedAny) {
		fputs(buf, outFile);
		// Update elapsed time and estimated time to completion
		float percentDone = (float)plussesPrinted / (float)totalPlusses;
		float seconds = (float) timer->Time();
		float estRemaining = seconds / percentDone - seconds;
		if (percentDone == 1.f)
			fprintf(outFile, " (%.1fs)       ", seconds);
		else
			fprintf(outFile, " (%.1fs|%.1fs)  ", seconds, max(0.f, estRemaining));
		fflush(outFile);
	}
}
void ProgressReporter::Done() const {
	while (plussesPrinted++ < totalPlusses)
		*curSpace++ = '+';
	fputs(buf, outFile);
	float seconds = (float) timer->Time();
	fprintf(outFile, " (%.1fs)       \n", seconds);
	fflush(outFile);
}
