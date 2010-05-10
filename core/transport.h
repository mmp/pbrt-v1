
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

#ifndef PBRT_TRANSPORT_H
#define PBRT_TRANSPORT_H
// transport.h*
#include "pbrt.h"
#include "primitive.h"
#include "color.h"
#include "light.h"
#include "reflection.h"
#include "sampling.h"
#include "material.h"
// Integrator Declarations
class COREDLL Integrator {
public:
	// Integrator Interface
	virtual ~Integrator();
	virtual Spectrum Li(const Scene *scene,
					    const RayDifferential &ray,
					    const Sample *sample,
					    float *alpha) const = 0;
	virtual void Preprocess(const Scene *scene) {
	}
	virtual void RequestSamples(Sample *sample,
	                            const Scene *scene) {
	}
};
class SurfaceIntegrator : public Integrator {
};
COREDLL Spectrum UniformSampleAllLights(const Scene *scene,
	const Point &p, const Normal &n, const Vector &wo,
	BSDF *bsdf, const Sample *sample,
	int *lightSampleOffset, int *bsdfSampleOffset,
	int *bsdfComponentOffset);
COREDLL Spectrum UniformSampleOneLight(const Scene *scene, const Point &p,
	const Normal &n, const Vector &wo, BSDF *bsdf,
	const Sample *sample, int lightSampleOffset = -1,
	int lightNumOffset = -1, int bsdfSampleOffset = -1,
	int bsdfComponentOffset = -1);
COREDLL Spectrum WeightedSampleOneLight(const Scene *scene, const Point &p,
	const Normal &n, const Vector &wo, BSDF *bsdf,
	const Sample *sample, int lightSampleOffset, int lightNumOffset,
	int bsdfSampleOffset, int bsdfComponentOffset, float *&avgY,
	float *&avgYsample, float *&cdf, float &overallAvgY);
#endif // PBRT_TRANSPORT_H
