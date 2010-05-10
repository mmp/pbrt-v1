
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

#ifndef PBRT_SAMPLING_H
#define PBRT_SAMPLING_H
// sampling.h*
#include "pbrt.h"
#include "geometry.h"
// Sampling Declarations
class COREDLL Sampler {
public:
	// Sampler Interface
	virtual ~Sampler();
	Sampler(int xstart, int xend,
	        int ystart, int yend,
			int spp);
	virtual bool GetNextSample(Sample *sample) = 0;
	int TotalSamples() const {
		return samplesPerPixel *
			(xPixelEnd - xPixelStart) *
			(yPixelEnd - yPixelStart);
	}
	virtual int RoundSize(int size) const = 0;
	// Sampler Public Data
	int xPixelStart, xPixelEnd, yPixelStart, yPixelEnd;
	int samplesPerPixel;
};
struct Sample {
	// Sample Public Methods
	Sample(SurfaceIntegrator *surf, VolumeIntegrator *vol,
		const Scene *scene);
	u_int Add1D(u_int num) {
		n1D.push_back(num);
		return n1D.size()-1;
	}
	u_int Add2D(u_int num) {
		n2D.push_back(num);
		return n2D.size()-1;
	}
	~Sample() {
		if (oneD != NULL) {
			FreeAligned(oneD[0]);
			FreeAligned(oneD);
		}
	}
	// Camera _Sample_ Data
	float imageX, imageY;
	float lensU, lensV;
	float time;
	// Integrator _Sample_ Data
	vector<u_int> n1D, n2D;
	float **oneD, **twoD;
};
COREDLL void StratifiedSample1D(float *samples,
					            int nsamples,
						        bool jitter = true);
COREDLL void StratifiedSample2D(float *samples,
                                int nx, int ny,
								bool jitter = true);
COREDLL void Shuffle(float *samp, int count, int dims);
COREDLL
void LatinHypercube(float *samples, int nSamples, int nDim);
inline double RadicalInverse(int n, int base) {
	double val = 0;
	double invBase = 1. / base, invBi = invBase;
	while (n > 0) {
		// Compute next digit of radical inverse
		int d_i = (n % base);
		val += d_i * invBi;
		n /= base;
		invBi *= invBase;
	}
	return val;
}
inline double FoldedRadicalInverse(int n, int base) {
	double val = 0;
	double invBase = 1.f/base, invBi = invBase;
	int modOffset = 0;
	while (val + base * invBi != val) {
		// Compute next digit of folded radical inverse
		int digit = ((n+modOffset) % base);
		val += digit * invBi;
		n /= base;
		invBi *= invBase;
		++modOffset;
	}
	return val;
}
inline float
	VanDerCorput(u_int n, u_int scramble = 0);
inline float
	Sobol2(u_int n, u_int scramble = 0);
inline float
	LarcherPillichshammer2(u_int n, u_int scramble = 0);
inline void
	Sample02(u_int n,
	            u_int scramble[2], float sample[2]);
class Filter {
public:
	// Filter Interface
	virtual ~Filter() { }
	Filter(float xw, float yw)
		: xWidth(xw), yWidth(yw), invXWidth(1.f/xw),
			invYWidth(1.f/yw) {
	}
	virtual float Evaluate(float x, float y) const = 0;
	// Filter Public Data
	const float xWidth, yWidth;
	const float invXWidth, invYWidth;
};
// Sampling Inline Functions
inline void Sample02(u_int n, u_int scramble[2],
                                  float sample[2]) {
	sample[0] = VanDerCorput(n, scramble[0]);
	sample[1] = Sobol2(n, scramble[1]);
}
inline float VanDerCorput(u_int n, u_int scramble) {
	n = (n << 16) | (n >> 16);
	n = ((n & 0x00ff00ff) << 8) | ((n & 0xff00ff00) >> 8);
	n = ((n & 0x0f0f0f0f) << 4) | ((n & 0xf0f0f0f0) >> 4);
	n = ((n & 0x33333333) << 2) | ((n & 0xcccccccc) >> 2);
	n = ((n & 0x55555555) << 1) | ((n & 0xaaaaaaaa) >> 1);
	n ^= scramble;
	return (float)n / (float)0x100000000LL;
}
inline float Sobol2(u_int n, u_int scramble) {
	for (u_int v = 1 << 31; n != 0; n >>= 1, v ^= v >> 1)
		if (n & 0x1) scramble ^= v;
	return (float)scramble / (float)0x100000000LL;
}
inline float
LarcherPillichshammer2(u_int n, u_int scramble) {
	for (u_int v = 1 << 31; n != 0; n >>= 1, v |= v >> 1)
		if (n & 0x1) scramble ^= v;
	return (float)scramble / (float)0x100000000LL;
}
inline void LDShuffleScrambled1D(int nSamples,
		int nPixel, float *samples) {
	u_int scramble = RandomUInt();
	for (int i = 0; i < nSamples * nPixel; ++i)
		samples[i] = VanDerCorput(i, scramble);
	for (int i = 0; i < nPixel; ++i)
		Shuffle(samples + i * nSamples, nSamples, 1);
	Shuffle(samples, nPixel, nSamples);
}
inline void LDShuffleScrambled2D(int nSamples,
		int nPixel, float *samples) {
	u_int scramble[2] = { RandomUInt(), RandomUInt() };
	for (int i = 0; i < nSamples * nPixel; ++i)
		Sample02(i, scramble, &samples[2*i]);
	for (int i = 0; i < nPixel; ++i)
		Shuffle(samples + 2 * i * nSamples, nSamples, 2);
	Shuffle(samples, nPixel, 2 * nSamples);
}
#endif // PBRT_SAMPLING_H
