
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

// bestcandidate.cpp*
#include "sampling.h"
#include "paramset.h"
#include "film.h"
// BestCandidate Sampling Constants
#define SQRT_SAMPLE_TABLE_SIZE 64
#define SAMPLE_TABLE_SIZE (SQRT_SAMPLE_TABLE_SIZE * \
                           SQRT_SAMPLE_TABLE_SIZE)
// BestCandidateSampler Declarations
class BestCandidateSampler : public Sampler {
public:
	// BestCandidateSampler Public Methods
	BestCandidateSampler(int xstart, int xend,
	                     int ystart, int yend,
						 int pixelsamples);
	~BestCandidateSampler() {
		delete[] strat2D;
		// so we leak on the individual elements of these arrays.  so it goes...
		delete[] oneDSamples;
		delete[] twoDSamples;
	}
	int RoundSize(int size) const {
		int root = Ceil2Int(sqrtf((float)size - .5f));
		return root*root;
	}
	bool GetNextSample(Sample *sample);
private:
	// BestCandidateSampler Private Data
	int tableOffset;
	float xTableCorner, yTableCorner, tableWidth;
	static const float sampleTable[SAMPLE_TABLE_SIZE][5];
	float **oneDSamples, **twoDSamples;
	int *strat2D;
	float sampleOffsets[3];
};
// BestCandidateSampler Method Definitions
BestCandidateSampler::
    BestCandidateSampler(int xstart, int xend,
		                 int ystart, int yend,
						 int pixelSamples)
	: Sampler(xstart, xend, ystart, yend, pixelSamples) {
	tableWidth =
		(float)SQRT_SAMPLE_TABLE_SIZE / sqrtf(pixelSamples);
	xTableCorner = float(xPixelStart) - tableWidth;
	yTableCorner = float(yPixelStart);
	tableOffset = SAMPLE_TABLE_SIZE;
	// _BestCandidateSampler_ constructor implementation
	oneDSamples = twoDSamples = NULL;
	strat2D = NULL;
}
#include "samplers/sampledata.cpp"
bool BestCandidateSampler::GetNextSample(Sample *sample) {
again:
	if (tableOffset == SAMPLE_TABLE_SIZE) {
		// Advance to next best-candidate sample table position
		tableOffset = 0;
		xTableCorner += tableWidth;
		if (xTableCorner >= xPixelEnd) {
			xTableCorner = float(xPixelStart);
			yTableCorner += tableWidth;
			if (yTableCorner >= yPixelEnd)
				return false;
		}
		if (!oneDSamples) {
			// Initialize sample tables and precompute _strat2D_ values
			oneDSamples = new float *[sample->n1D.size()];
			for (u_int i = 0; i < sample->n1D.size(); ++i) {
				oneDSamples[i] = (sample->n1D[i] == 1) ?
					new float[SAMPLE_TABLE_SIZE] : NULL;
			}
			twoDSamples = new float *[sample->n2D.size()];
			strat2D = new int[sample->n2D.size()];
			for (u_int i = 0; i < sample->n2D.size(); ++i) {
				twoDSamples[i] = (sample->n2D[i] == 1) ?
					new float[2 * SAMPLE_TABLE_SIZE] : NULL;
				strat2D[i] =
					Ceil2Int(sqrtf((float)sample->n2D[i] - .5f));
			}
		}
		// Update sample shifts
		for (int i = 0; i < 3; ++i)
			sampleOffsets[i] = RandomFloat();
		// Generate _SAMPLE\_TABLE\_SIZE_-sized tables for single samples
		for (u_int i = 0; i < sample->n1D.size(); ++i)
			if (sample->n1D[i] == 1)
				LDShuffleScrambled1D(SAMPLE_TABLE_SIZE, 1,
				                     oneDSamples[i]);
		for (u_int i = 0; i < sample->n2D.size(); ++i)
			if (sample->n2D[i] == 1)
				LDShuffleScrambled2D(SAMPLE_TABLE_SIZE, 1,
				                     twoDSamples[i]);
	}
	// Compute raster sample from table
	#define WRAP(x) ((x) > 1 ? ((x)-1) : (x))
	sample->imageX = xTableCorner + tableWidth *
		sampleTable[tableOffset][0];
	sample->imageY = yTableCorner + tableWidth *
		sampleTable[tableOffset][1];
	sample->time  = WRAP(sampleOffsets[0] +
		sampleTable[tableOffset][2]);
	sample->lensU = WRAP(sampleOffsets[1] +
		sampleTable[tableOffset][3]);
	sample->lensV = WRAP(sampleOffsets[2] +
		sampleTable[tableOffset][4]);
	// Check sample against crop window, goto _again_ if outside
	if (sample->imageX <  xPixelStart ||
	    sample->imageX >= xPixelEnd   ||
	    sample->imageY <  yPixelStart ||
	    sample->imageY >= yPixelEnd) {
		++tableOffset;
		goto again;
	}
	// Compute integrator samples for best-candidate sample
	for (u_int i = 0; i < sample->n1D.size(); ++i) {
		if (sample->n1D[i] == 1)
			sample->oneD[i][0] = oneDSamples[i][tableOffset];
		else
			StratifiedSample1D(sample->oneD[i], sample->n1D[i]);
	}
	for (u_int i = 0; i < sample->n2D.size(); ++i) {
		if (sample->n2D[i] == 1) {
		   sample->twoD[i][0] = twoDSamples[i][2*tableOffset];
		   sample->twoD[i][1] = twoDSamples[i][2*tableOffset+1];
		}
		else {
			StratifiedSample2D(sample->twoD[i],
			                   strat2D[i],
							   strat2D[i]);
		}
	}
	++tableOffset;
	return true;
}
extern "C" DLLEXPORT Sampler *CreateSampler(const ParamSet &params, const Film *film) {
	// Initialize common sampler parameters
	int xstart, xend, ystart, yend;
	film->GetSampleExtent(&xstart, &xend, &ystart, &yend);
	int nsamp = params.FindOneInt("pixelsamples", 4);
	return new BestCandidateSampler(xstart, xend, ystart, yend, nsamp);
}
