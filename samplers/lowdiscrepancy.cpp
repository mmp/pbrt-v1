
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

// lowdiscrepancy.cpp*
#include "sampling.h"
#include "paramset.h"
#include "film.h"
// LDSampler Declarations
class LDSampler : public Sampler {
public:
	// LDSampler Public Methods
	LDSampler(int xstart, int xend,
	          int ystart, int yend,
			  int nsamp);
	~LDSampler() {
		delete[] imageSamples;
		for (int i = 0; i < n1D; ++i)
			delete[] oneDSamples[i];
		for (int i = 0; i < n2D; ++i)
			delete[] twoDSamples[i];
		delete[] oneDSamples;
		delete[] twoDSamples;
	}
	int RoundSize(int size) const {
		return RoundUpPow2(size);
	}
	bool GetNextSample(Sample *sample);
private:
	// LDSampler Private Data
	int xPos, yPos, pixelSamples;
	int samplePos;
	float *imageSamples, *lensSamples, *timeSamples;
	float **oneDSamples, **twoDSamples;
	int n1D, n2D;
};
// LDSampler Method Definitions
LDSampler::LDSampler(int xstart, int xend,
		int ystart, int yend, int ps)
	: Sampler(xstart, xend, ystart, yend, RoundUpPow2(ps)) {
	xPos = xPixelStart - 1;
	yPos = yPixelStart;
	if (!IsPowerOf2(ps)) {
		Warning("Pixel samples being"
		        " rounded up to power of 2");
		pixelSamples = RoundUpPow2(ps);
	}
	else
		pixelSamples = ps;
	samplePos = pixelSamples;
	oneDSamples = twoDSamples = NULL;
	imageSamples = new float[5*pixelSamples];
	lensSamples = imageSamples + 2*pixelSamples;
	timeSamples = imageSamples + 4*pixelSamples;
	n1D = n2D = 0;
}
bool LDSampler::GetNextSample(Sample *sample) {
	if (!oneDSamples) {
		// Allocate space for pixel's low-discrepancy sample tables
		oneDSamples = new float *[sample->n1D.size()];
		n1D = sample->n1D.size();
		for (u_int i = 0; i < sample->n1D.size(); ++i)
			oneDSamples[i] = new float[sample->n1D[i] *
		                               pixelSamples];
		twoDSamples = new float *[sample->n2D.size()];
		n2D = sample->n2D.size();
		for (u_int i = 0; i < sample->n2D.size(); ++i)
			twoDSamples[i] = new float[2 * sample->n2D[i] *
		                               pixelSamples];
	}
	if (samplePos == pixelSamples) {
		// Advance to next pixel for low-discrepancy sampling
		if (++xPos == xPixelEnd) {
			xPos = xPixelStart;
			++yPos;
		}
		if (yPos == yPixelEnd)
			return false;
		samplePos = 0;
		// Generate low-discrepancy samples for pixel
		LDShuffleScrambled2D(1, pixelSamples, imageSamples);
		LDShuffleScrambled2D(1, pixelSamples, lensSamples);
		LDShuffleScrambled1D(1, pixelSamples, timeSamples);
		for (u_int i = 0; i < sample->n1D.size(); ++i)
			LDShuffleScrambled1D(sample->n1D[i], pixelSamples,
				oneDSamples[i]);
		for (u_int i = 0; i < sample->n2D.size(); ++i)
			LDShuffleScrambled2D(sample->n2D[i], pixelSamples,
				twoDSamples[i]);
	}
	// Copy low-discrepancy samples from tables
	sample->imageX = xPos + imageSamples[2*samplePos];
	sample->imageY = yPos + imageSamples[2*samplePos+1];
	sample->time = timeSamples[samplePos];
	sample->lensU = lensSamples[2*samplePos];
	sample->lensV = lensSamples[2*samplePos+1];
	for (u_int i = 0; i < sample->n1D.size(); ++i) {
		int startSamp = sample->n1D[i] * samplePos;
		for (u_int j = 0; j < sample->n1D[i]; ++j)
			sample->oneD[i][j] = oneDSamples[i][startSamp+j];
	}
	for (u_int i = 0; i < sample->n2D.size(); ++i) {
		int startSamp = 2 * sample->n2D[i] * samplePos;
		for (u_int j = 0; j < 2*sample->n2D[i]; ++j)
			sample->twoD[i][j] = twoDSamples[i][startSamp+j];
	}
	++samplePos;
	return true;
}
extern "C" DLLEXPORT Sampler *CreateSampler(const ParamSet &params, const Film *film) {
	// Initialize common sampler parameters
	int xstart, xend, ystart, yend;
	film->GetSampleExtent(&xstart, &xend, &ystart, &yend);
	int nsamp = params.FindOneInt("pixelsamples", 4);
	return new LDSampler(xstart, xend, ystart, yend, nsamp);
}
