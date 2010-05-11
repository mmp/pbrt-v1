
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
