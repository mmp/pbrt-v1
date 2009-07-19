
/*
 * pbrt source code Copyright(c) 1998-2005 Matt Pharr and Greg Humphreys
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 * (See file License.txt for complete license)
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
