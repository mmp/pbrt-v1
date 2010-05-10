
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

// directlighting.cpp*
#include "pbrt.h"
#include "transport.h"
#include "scene.h"
// DirectLighting Declarations
enum LightStrategy { SAMPLE_ALL_UNIFORM, SAMPLE_ONE_UNIFORM,
	SAMPLE_ONE_WEIGHTED  // NOBOOK
};
class DirectLighting : public SurfaceIntegrator {
public:
	// DirectLighting Public Methods
	DirectLighting(LightStrategy ls, int md);
	~DirectLighting();
	Spectrum Li(const Scene *scene, const RayDifferential &ray, const Sample *sample,
		float *alpha) const;
	void RequestSamples(Sample *sample, const Scene *scene) {
		if (strategy == SAMPLE_ALL_UNIFORM) {
			// Allocate and request samples for sampling all lights
			u_int nLights = scene->lights.size();
			lightSampleOffset = new int[nLights];
			bsdfSampleOffset = new int[nLights];
			bsdfComponentOffset = new int[nLights];
			for (u_int i = 0; i < nLights; ++i) {
				const Light *light = scene->lights[i];
				int lightSamples =
					scene->sampler->RoundSize(light->nSamples);
				lightSampleOffset[i] = sample->Add2D(lightSamples);
				bsdfSampleOffset[i] = sample->Add2D(lightSamples);
				bsdfComponentOffset[i] = sample->Add1D(lightSamples);
			}
			lightNumOffset = -1;
		}
		else {
			// Allocate and request samples for sampling one light
			lightSampleOffset = new int[1];
			lightSampleOffset[0] = sample->Add2D(1);
			lightNumOffset = sample->Add1D(1);
			bsdfSampleOffset = new int[1];
			bsdfSampleOffset[0] = sample->Add2D(1);
			bsdfComponentOffset = new int[1];
			bsdfComponentOffset[0] = sample->Add1D(1);
		}
	}
private:
	// DirectLighting Private Data
	LightStrategy strategy;
	mutable int rayDepth; // NOBOOK
	int maxDepth; // NOBOOK
	// Declare sample parameters for light source sampling
	int *lightSampleOffset, lightNumOffset;
	int *bsdfSampleOffset, *bsdfComponentOffset;
	mutable float *avgY, *avgYsample, *cdf;
	mutable float overallAvgY;
};
// DirectLighting Method Definitions
DirectLighting::~DirectLighting() {
	delete[] avgY;
	delete[] avgYsample;
	delete[] cdf;
}
DirectLighting::DirectLighting(LightStrategy st, int md) {
	maxDepth = md;
	rayDepth = 0;
	strategy = st;
	avgY = avgYsample = cdf = NULL;
	overallAvgY = 0.;
}
Spectrum DirectLighting::Li(const Scene *scene,
		const RayDifferential &ray, const Sample *sample,
		float *alpha) const {
	Intersection isect;
	Spectrum L(0.);
	if (scene->Intersect(ray, &isect)) {
		if (alpha) *alpha = 1.;
		// Evaluate BSDF at hit point
		BSDF *bsdf = isect.GetBSDF(ray);
		Vector wo = -ray.d;
		const Point &p = bsdf->dgShading.p;
		const Normal &n = bsdf->dgShading.nn;
		// Compute emitted light if ray hit an area light source
		L += isect.Le(wo);
		// Compute direct lighting for _DirectLighting_ integrator
		if (scene->lights.size() > 0) {
			// Apply direct lighting strategy
			switch (strategy) {
				case SAMPLE_ALL_UNIFORM:
					L += UniformSampleAllLights(scene, p, n, wo, bsdf,
						sample, lightSampleOffset, bsdfSampleOffset,
						bsdfComponentOffset);
					break;
				case SAMPLE_ONE_UNIFORM:
					L += UniformSampleOneLight(scene, p, n, wo, bsdf,
						sample, lightSampleOffset[0], lightNumOffset,
						bsdfSampleOffset[0], bsdfComponentOffset[0]);
					break;
				case SAMPLE_ONE_WEIGHTED:
					L += WeightedSampleOneLight(scene, p, n, wo, bsdf,
						sample, lightSampleOffset[0], lightNumOffset,
						bsdfSampleOffset[0], bsdfComponentOffset[0], avgY,
						avgYsample, cdf, overallAvgY);
					break;
			}
		}
		if (rayDepth++ < maxDepth) {
			Vector wi;
			// Trace rays for specular reflection and refraction
			Spectrum f = bsdf->Sample_f(wo, &wi,
				BxDFType(BSDF_REFLECTION | BSDF_SPECULAR));
			if (!f.Black() && AbsDot(wi, n) > 0.f) {
				// Compute ray differential _rd_ for specular reflection
				RayDifferential rd(p, wi);
				rd.hasDifferentials = true;
				rd.rx.o = p + isect.dg.dpdx;
				rd.ry.o = p + isect.dg.dpdy;
				// Compute differential reflected directions
				Normal dndx = bsdf->dgShading.dndu * bsdf->dgShading.dudx +
					bsdf->dgShading.dndv * bsdf->dgShading.dvdx;
				Normal dndy = bsdf->dgShading.dndu * bsdf->dgShading.dudy +
					bsdf->dgShading.dndv * bsdf->dgShading.dvdy;
				Vector dwodx = -ray.rx.d - wo, dwody = -ray.ry.d - wo;
				float dDNdx = Dot(dwodx, n) + Dot(wo, dndx);
				float dDNdy = Dot(dwody, n) + Dot(wo, dndy);
				rd.rx.d = wi -
				          dwodx + 2 * Vector(Dot(wo, n) * dndx +
						  dDNdx * n);
				rd.ry.d = wi -
				          dwody + 2 * Vector(Dot(wo, n) * dndy +
						  dDNdy * n);
				L += scene->Li(rd, sample) * f * AbsDot(wi, n);
			}
			f = bsdf->Sample_f(wo, &wi,
				BxDFType(BSDF_TRANSMISSION | BSDF_SPECULAR));
			if (!f.Black() && AbsDot(wi, n) > 0.f) {
				// Compute ray differential _rd_ for specular transmission
				RayDifferential rd(p, wi);
				rd.hasDifferentials = true;
				rd.rx.o = p + isect.dg.dpdx;
				rd.ry.o = p + isect.dg.dpdy;
				
				float eta = bsdf->eta;
				Vector w = -wo;
				if (Dot(wo, n) < 0) eta = 1.f / eta;
				
				Normal dndx = bsdf->dgShading.dndu * bsdf->dgShading.dudx + bsdf->dgShading.dndv * bsdf->dgShading.dvdx;
				Normal dndy = bsdf->dgShading.dndu * bsdf->dgShading.dudy + bsdf->dgShading.dndv * bsdf->dgShading.dvdy;
				
				Vector dwodx = -ray.rx.d - wo, dwody = -ray.ry.d - wo;
				float dDNdx = Dot(dwodx, n) + Dot(wo, dndx);
				float dDNdy = Dot(dwody, n) + Dot(wo, dndy);
				
				float mu = eta * Dot(w, n) - Dot(wi, n);
				float dmudx = (eta - (eta*eta*Dot(w,n))/Dot(wi, n)) * dDNdx;
				float dmudy = (eta - (eta*eta*Dot(w,n))/Dot(wi, n)) * dDNdy;
				
				rd.rx.d = wi + eta * dwodx - Vector(mu * dndx + dmudx * n);
				rd.ry.d = wi + eta * dwody - Vector(mu * dndy + dmudy * n);
				L += scene->Li(rd, sample) * f * AbsDot(wi, n);
			}
		}
		--rayDepth;
	}
	else {
		// Handle ray with no intersection
		if (alpha) *alpha = 0.;
		for (u_int i = 0; i < scene->lights.size(); ++i)
			L += scene->lights[i]->Le(ray);
		if (alpha && !L.Black()) *alpha = 1.;
		return L;
	}
	return L;
}
extern "C" DLLEXPORT SurfaceIntegrator *CreateSurfaceIntegrator(const ParamSet &params) {
	int maxDepth = params.FindOneInt("maxdepth", 5);
	LightStrategy strategy;
	string st = params.FindOneString("strategy", "all");
	if (st == "one") strategy = SAMPLE_ONE_UNIFORM;
	else if (st == "all") strategy = SAMPLE_ALL_UNIFORM;
	else if (st == "weighted") strategy = SAMPLE_ONE_WEIGHTED;
	else {
		Warning("Strategy \"%s\" for direct lighting unknown. "
			"Using \"all\".", st.c_str());
		strategy = SAMPLE_ALL_UNIFORM;
	}
	return new DirectLighting(strategy, maxDepth);
}
