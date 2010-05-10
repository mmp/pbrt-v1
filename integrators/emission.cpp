
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

// emission.cpp*
#include "volume.h"
#include "transport.h"
#include "scene.h"
// EmissionIntegrator Declarations
class EmissionIntegrator : public VolumeIntegrator {
public:
	// EmissionIntegrator Public Methods
	EmissionIntegrator(float ss) { stepSize = ss; }
	void RequestSamples(Sample *sample, const Scene *scene);
	Spectrum Transmittance(const Scene *, const Ray &ray, const Sample *sample, float *alpha) const;
	Spectrum Li(const Scene *, const RayDifferential &ray, const Sample *sample, float *alpha) const;
private:
	// EmissionIntegrator Private Data
	float stepSize;
	int tauSampleOffset, scatterSampleOffset;
};
// EmissionIntegrator Method Definitions
void EmissionIntegrator::RequestSamples(Sample *sample,
		const Scene *scene) {
	tauSampleOffset = sample->Add1D(1);
	scatterSampleOffset = sample->Add1D(1);
}
Spectrum
	EmissionIntegrator::Transmittance(const Scene *scene,
		const Ray &ray, const Sample *sample,
		float *alpha) const {
	if (!scene->volumeRegion) return Spectrum(1.f);
	float step = sample ? stepSize : 4.f * stepSize;
	float offset =
		sample ? sample->oneD[tauSampleOffset][0] :
		RandomFloat();
	Spectrum tau =
		scene->volumeRegion->Tau(ray, step, offset);
	return Exp(-tau);
}
Spectrum EmissionIntegrator::Li(const Scene *scene,
		const RayDifferential &ray, const Sample *sample,
		float *alpha) const {
	VolumeRegion *vr = scene->volumeRegion;
	float t0, t1;
	if (!vr || !vr->IntersectP(ray, &t0, &t1) || (t1-t0) == 0.f) return 0.f;
	// Do emission-only volume integration in _vr_
	Spectrum Lv(0.);
	// Prepare for volume integration stepping
	int N = Ceil2Int((t1-t0) / stepSize);
	float step = (t1 - t0) / N;
	Spectrum Tr(1.f);
	Point p = ray(t0), pPrev;
	Vector w = -ray.d;
	if (sample)
		t0 += sample->oneD[scatterSampleOffset][0] * step;
	else
		t0 += RandomFloat() * step;
	for (int i = 0; i < N; ++i, t0 += step) {
		// Advance to sample at _t0_ and update _T_
		pPrev = p;
		p = ray(t0);
		Spectrum stepTau = vr->Tau(Ray(pPrev, p - pPrev, 0, 1),
			.5f * stepSize, RandomFloat());
		Tr *= Exp(-stepTau);
		// Possibly terminate raymarching if transmittance is small
		if (Tr.y() < 1e-3) {
			const float continueProb = .5f;
			if (RandomFloat() > continueProb) break;
			Tr /= continueProb;
		}
		// Compute emission-only source term at _p_
		Lv += Tr * vr->Lve(p, w);
	}
	return Lv * step;
}
extern "C" DLLEXPORT VolumeIntegrator *CreateVolumeIntegrator(const ParamSet &params) {
	float stepSize  = params.FindOneFloat("stepsize", 1.f);
	return new EmissionIntegrator(stepSize);
}
