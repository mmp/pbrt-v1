
/*
 * pbrt source code Copyright(c) 1998-2007 Matt Pharr and Greg Humphreys
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 * (See file License.txt for complete license)
 */

// sampling.cpp*
#include "pbrt.h"
#include "sampling.h"
#include "transport.h"
#include "volume.h"
// Sampler Method Definitions
Sampler::~Sampler() {
}
Sampler::Sampler(int xstart, int xend, int ystart, int yend,
		int spp) {
	xPixelStart = xstart;
	xPixelEnd = xend;
	yPixelStart = ystart;
	yPixelEnd = yend;
	samplesPerPixel = spp;
}
// Sample Method Definitions
Sample::Sample(SurfaceIntegrator *surf,
		VolumeIntegrator *vol, const Scene *scene) {
	surf->RequestSamples(this, scene);
	vol->RequestSamples(this, scene);
	// Allocate storage for sample pointers
	int nPtrs = n1D.size() + n2D.size();
	if (!nPtrs) {
		oneD = twoD = NULL;
		return;
	}
	oneD = (float **)AllocAligned(nPtrs * sizeof(float *));
	twoD = oneD + n1D.size();
	// Compute total number of sample values needed
	int totSamples = 0;
	for (u_int i = 0; i < n1D.size(); ++i)
		totSamples += n1D[i];
	for (u_int i = 0; i < n2D.size(); ++i)
		totSamples += 2 * n2D[i];
	// Allocate storage for sample values
	float *mem = (float *)AllocAligned(totSamples *
		sizeof(float));
	for (u_int i = 0; i < n1D.size(); ++i) {
		oneD[i] = mem;
		mem += n1D[i];
	}
	for (u_int i = 0; i < n2D.size(); ++i) {
		twoD[i] = mem;
		mem += 2 * n2D[i];
	}
}
// Sampling Function Definitions
COREDLL void StratifiedSample1D(float *samp, int nSamples,
		bool jitter) {
	float invTot = 1.f / nSamples;
	for (int i = 0;  i < nSamples; ++i) {
		float delta = jitter ? RandomFloat() : 0.5f;
		*samp++ = (i + delta) * invTot;
	}
}
COREDLL void StratifiedSample2D(float *samp, int nx, int ny,
		bool jitter) {
	float dx = 1.f / nx, dy = 1.f / ny;
	for (int y = 0; y < ny; ++y)
		for (int x = 0; x < nx; ++x) {
			float jx = jitter ? RandomFloat() : 0.5f;
			float jy = jitter ? RandomFloat() : 0.5f;
			*samp++ = (x + jx) * dx;
			*samp++ = (y + jy) * dy;
		}
}
COREDLL void Shuffle(float *samp, int count, int dims) {
	for (int i = 0; i < count; ++i) {
		u_int other = RandomUInt() % count;
		for (int j = 0; j < dims; ++j)
			swap(samp[dims*i + j], samp[dims*other + j]);
	}
}
COREDLL void LatinHypercube(float *samples,
                             int nSamples, int nDim) {
	// Generate LHS samples along diagonal
	float delta = 1.f / nSamples;
	for (int i = 0; i < nSamples; ++i)
		for (int j = 0; j < nDim; ++j)
			samples[nDim * i + j] = (i + RandomFloat()) * delta;
	// Permute LHS samples in each dimension
	for (int i = 0; i < nDim; ++i) {
		for (int j = 0; j < nSamples; ++j) {
			u_int other = RandomUInt() % nSamples;
			swap(samples[nDim * j + i],
			     samples[nDim * other + i]);
		}
	}
}
