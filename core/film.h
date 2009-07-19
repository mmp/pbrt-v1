
/*
 * pbrt source code Copyright(c) 1998-2007 Matt Pharr and Greg Humphreys
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 * (See file License.txt for complete license)
 */

#ifndef PBRT_FILM_H
#define PBRT_FILM_H
// film.h*
#include "pbrt.h"
// Film Declarations
class Film {
public:
	// Film Interface
	Film(int xres, int yres)
		: xResolution(xres), yResolution(yres) {
	}
	virtual ~Film() { }
	virtual void AddSample(const Sample &sample, const Ray &ray,
		const Spectrum &L, float alpha) = 0;
	virtual void WriteImage() = 0;
	virtual void GetSampleExtent(int *xstart,
		int *xend, int *ystart, int *yend) const = 0;
	// Film Public Data
	const int xResolution, yResolution;
};
// Image Pipeline Declarations
extern COREDLL void ApplyImagingPipeline(float *rgb,
	int xResolution, int yResolution,
	float *yWeight = NULL,
	float bloomRadius = .2f, float bloomWeight = 0.f,
	const char *tonemap = NULL,
	const ParamSet *toneMapParams = NULL,
	float gamma = 2.2, float dither = 0.5f,
	int maxDisplayValue = 255);
#endif // PBRT_FILM_H
