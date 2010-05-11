
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

// random.cpp*
#include "sampling.h"
#include "paramset.h"
#include "film.h"
class RandomSampler : public Sampler {
public:
	RandomSampler(int xstart, int xend, int ystart,
		int yend, int xs, int ys);
	~RandomSampler() {
		FreeAligned(imageSamples);
	}
	bool GetNextSample(Sample *sample);
	int RoundSize(int sz) const { return sz; }
private:
	// RandomSampler Private Data
	bool jitterSamples;
	int xPos, yPos, xPixelSamples, yPixelSamples;
	float *imageSamples, *lensSamples, *timeSamples;
	int samplePos;
};

RandomSampler::RandomSampler(int xstart, int xend,
		int ystart, int yend, int xs, int ys)
	: Sampler(xstart, xend, ystart, yend, xs * ys) {
	xPos = xPixelStart;
	yPos = yPixelStart;
	xPixelSamples = xs;
	yPixelSamples = ys;
	// Get storage for a pixel's worth of stratified samples
	imageSamples = (float *)AllocAligned(5 * xPixelSamples *
		yPixelSamples * sizeof(float));
	lensSamples = imageSamples +
	              2 * xPixelSamples * yPixelSamples;
	timeSamples = lensSamples +
	              2 * xPixelSamples * yPixelSamples;

	for (int i = 0;
	     i < 5 * xPixelSamples * yPixelSamples;
		 ++i) {
		imageSamples[i] = RandomFloat();
	}

	// Shift image samples to pixel coordinates
	for (int o = 0;
	     o < 2 * xPixelSamples * yPixelSamples;
		 o += 2) {
		imageSamples[o]   += xPos;
		imageSamples[o+1] += yPos;
	}
	samplePos = 0;
}

bool RandomSampler::GetNextSample(Sample *sample) {
	// Compute new set of samples if needed for next pixel
	if (samplePos == xPixelSamples * yPixelSamples) {
		// Advance to next pixel for stratified sampling
		if (++xPos == xPixelEnd) {
			xPos = xPixelStart;
			++yPos;
		}
		if (yPos == yPixelEnd)
			return false;

		for (int i = 0;
		     i < 5 * xPixelSamples * yPixelSamples;
			 ++i) {
			imageSamples[i] = RandomFloat();
		}

		// Shift image samples to pixel coordinates
		for (int o = 0;
		     o < 2 * xPixelSamples * yPixelSamples;
			 o += 2) {
			imageSamples[o]   += xPos;
			imageSamples[o+1] += yPos;
		}
		samplePos = 0;
	}
	// Return next \mono{RandomSampler} sample point
	sample->imageX = imageSamples[2*samplePos];
	sample->imageY = imageSamples[2*samplePos+1];
	sample->lensU = lensSamples[2*samplePos];
	sample->lensV = lensSamples[2*samplePos+1];
	sample->time = timeSamples[samplePos];
	// Generate stratified samples for integrators
	for (u_int i = 0; i < sample->n1D.size(); ++i)
		for (u_int j = 0; j < sample->n1D[i]; ++j)
			sample->oneD[i][j] = RandomFloat();
	for (u_int i = 0; i < sample->n2D.size(); ++i)
		for (u_int j = 0; j < 2*sample->n2D[i]; ++j)
			sample->twoD[i][j] = RandomFloat();
	++samplePos;
	return true;
}

extern "C" DLLEXPORT
Sampler *CreateSampler(const ParamSet &params,
                       const Film *film) {
	int xsamp = params.FindOneInt("xsamples", 2);
	int ysamp = params.FindOneInt("ysamples", 2);
	int xstart, xend, ystart, yend;
	film->GetSampleExtent(&xstart, &xend, &ystart, &yend);
	return new RandomSampler(xstart, xend,
	                         ystart, yend,
							 xsamp, ysamp);
}
