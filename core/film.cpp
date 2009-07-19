
/*
 * pbrt source code Copyright(c) 1998-2007 Matt Pharr and Greg Humphreys
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 * (See file License.txt for complete license)
 */

// film.cpp*
#include "film.h"
#include "dynload.h"
#include "paramset.h"
#include "tonemap.h"
// Image Pipeline Function Definitions
void ApplyImagingPipeline(float *rgb, int xResolution,
		int yResolution, float *yWeight,
		float bloomRadius, float bloomWeight,
		const char *toneMapName,
		const ParamSet *toneMapParams,
		float gamma, float dither, int maxDisplayValue) {
	int nPix = xResolution * yResolution ;
	// Possibly apply bloom effect to image
	if (bloomRadius > 0.f && bloomWeight > 0.f) {
		// Compute image-space extent of bloom effect
		int bloomSupport = Float2Int(bloomRadius *
			max(xResolution, yResolution));
		int bloomWidth = bloomSupport / 2;
		// Initialize bloom filter table
		float *bloomFilter = new float[bloomWidth * bloomWidth];
		for (int i = 0; i < bloomWidth * bloomWidth; ++i) {
			float dist = sqrtf(float(i)) / float(bloomWidth);
			bloomFilter[i] = powf(max(0.f, 1.f - dist), 4.f);
		}
		// Apply bloom filter to image pixels
		float *bloomImage = new float[3*nPix];
		ProgressReporter prog(yResolution, "Bloom filter"); //NOBOOK
		for (int y = 0; y < yResolution; ++y) {
			for (int x = 0; x < xResolution; ++x) {
				// Compute bloom for pixel _(x,y)_
				// Compute extent of pixels contributing bloom
				int x0 = max(0, x - bloomWidth);
				int x1 = min(x + bloomWidth, xResolution - 1);
				int y0 = max(0, y - bloomWidth);
				int y1 = min(y + bloomWidth, yResolution - 1);
				int offset = y * xResolution + x;
				float sumWt = 0.;
				for (int by = y0; by <= y1; ++by)
					for (int bx = x0; bx <= x1; ++bx) {
						// Accumulate bloom from pixel $(bx,by)$
						int dx = x - bx, dy = y - by;
						if (dx == 0 && dy == 0) continue;
						int dist2 = dx*dx + dy*dy;
						if (dist2 < bloomWidth * bloomWidth) {
							int bloomOffset = bx + by * xResolution;
							float wt = bloomFilter[dist2];
							sumWt += wt;
							for (int j = 0; j < 3; ++j)
								bloomImage[3*offset+j] += wt * rgb[3*bloomOffset+j];
						}
					}
				bloomImage[3*offset  ] /= sumWt;
				bloomImage[3*offset+1] /= sumWt;
				bloomImage[3*offset+2] /= sumWt;
			}
			prog.Update(); //NOBOOK
		}
		prog.Done(); //NOBOOK
		// Mix bloom effect into each pixel
		for (int i = 0; i < 3 * nPix; ++i)
			rgb[i] = Lerp(bloomWeight, rgb[i], bloomImage[i]);
		// Free memory allocated for bloom effect
		delete[] bloomFilter;
		delete[] bloomImage;
	}
	// Apply tone reproduction to image
	ToneMap *toneMap = NULL;
	if (toneMapName)
		toneMap = MakeToneMap(toneMapName,
			toneMapParams ? *toneMapParams : ParamSet());
	if (toneMap) {
		float maxDisplayY = 100.f;
		float *scale = new float[nPix], *lum = new float[nPix];
		// Compute pixel luminance values
		float stdYWeight[3] = { 0.212671f, 0.715160f, 0.072169f };
		if (!yWeight) yWeight = stdYWeight;
		for (int i = 0; i < nPix; ++i)
			lum[i] = 683.f * (yWeight[0] * rgb[3*i] +
				yWeight[1] * rgb[3*i+1] + yWeight[2] * rgb[3*i+2]);
		toneMap->Map(lum, xResolution, yResolution,
			maxDisplayY, scale);
		// Apple scale to pixels for tone mapping and map to $[0,1]$
		float displayTo01 = 683.f / maxDisplayY;
		for (int i = 0; i < xResolution * yResolution; ++i) {
			rgb[3*i  ] *= scale[i] * displayTo01;
			rgb[3*i+1] *= scale[i] * displayTo01;
			rgb[3*i+2] *= scale[i] * displayTo01;
		}
		delete[] scale;
		delete[] lum;
	}
	// Handle out-of-gamut RGB values
	for (int i = 0; i < nPix; ++i) {
		float m = max(rgb[3*i], max(rgb[3*i+1], rgb[3*i+2]));
		if (m > 1.f)
			for (int j = 0; j < 3; ++j)
				rgb[3*i+j] /= m;
	}
	// Apply gamma correction to image
	if (gamma != 1.f) {
		float invGamma = 1.f / gamma;
		for (int i = 0; i < 3*nPix; ++i)
			rgb[i] = powf(rgb[i], invGamma);
	}
	// Map image to display range
	for (int i = 0; i < 3*nPix; ++i)
		rgb[i] *= maxDisplayValue;
	// Dither image
	if (dither > 0.f)
		for (int i = 0; i < 3*nPix; ++i)
			rgb[i] += 2.f * dither * (RandomFloat() - .5f);
}
