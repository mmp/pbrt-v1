
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

// igi.cpp*
#include "pbrt.h"
#include "transport.h"
#include "scene.h"
#include "mc.h"
#include "sampling.h"

// IGI Local Structures
struct VirtualLight {
	VirtualLight() { }
	VirtualLight(const Point &pp, const Normal &nn, const Spectrum &le)
		: p(pp), n(nn), Le(le) { }
	Point p;
	Normal n;
	Spectrum Le;
};

class IGIIntegrator : public SurfaceIntegrator {
public:
	// IGIIntegrator Public Methods
	IGIIntegrator(int nl, int ns, float md, float rrt, float is);
	Spectrum Li(const Scene *scene, const RayDifferential &ray,
		const Sample *sample, float *alpha) const;
	void RequestSamples(Sample *sample, const Scene *scene);
	void Preprocess(const Scene *);

private:
	// IGI Private Data
	u_int nLightPaths, nLightSets;
	vector<VirtualLight> *virtualLights;
	mutable int specularDepth;
	int maxSpecularDepth;
	float minDist2, rrThreshold, indirectScale;
	int vlSetOffset;

	int *lightSampleOffset, lightNumOffset;
	int *bsdfSampleOffset, *bsdfComponentOffset;
};

// IGIIntegrator Implementation
IGIIntegrator::IGIIntegrator(int nl, int ns, float md,
		float rrt, float is) {
	nLightPaths = RoundUpPow2((u_int)nl);
	nLightSets = RoundUpPow2((u_int)ns);
	minDist2 = md * md;
	rrThreshold = rrt;
	indirectScale = is;
	maxSpecularDepth = 5;
	specularDepth = 0;
	virtualLights = new vector<VirtualLight>[nLightSets];
}
void IGIIntegrator::RequestSamples(Sample *sample,
                                   const Scene *scene) {
	// Request samples for area light sampling
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
	vlSetOffset = sample->Add1D(1);
}
void IGIIntegrator::Preprocess(const Scene *scene) {
	if (scene->lights.size() == 0) return;
	// Compute samples for emitted rays from lights
	float *lightNum = new float[nLightPaths * nLightSets];
	float *lightSamp0 = new float[2 * nLightPaths *	nLightSets];
	float *lightSamp1 = new float[2 * nLightPaths * nLightSets];
	LDShuffleScrambled1D(nLightPaths, nLightSets, lightNum);
	LDShuffleScrambled2D(nLightPaths, nLightSets, lightSamp0);
	LDShuffleScrambled2D(nLightPaths, nLightSets, lightSamp1);
	// Precompute information for light sampling densities
	int nLights = int(scene->lights.size());
	float *lightPower = (float *)alloca(nLights * sizeof(float));
	float *lightCDF = (float *)alloca((nLights+1) * sizeof(float));
	for (int i = 0; i < nLights; ++i)
		lightPower[i] = scene->lights[i]->Power(scene).y();
	float totalPower;
	ComputeStep1dCDF(lightPower, nLights, &totalPower, lightCDF);
	for (u_int s = 0; s < nLightSets; ++s) {
		for (u_int i = 0; i < nLightPaths; ++i) {
			// Follow path _i_ from light to create virtual lights
			int sampOffset = s*nLightPaths + i;
			// Choose light source to trace path from
			float lightPdf;
			int lNum = Floor2Int(SampleStep1d(lightPower, lightCDF,
				totalPower, nLights, lightNum[sampOffset], &lightPdf) * nLights);
//			fprintf(stderr, "samp %f -> num %d\n", lightNum[sampOffset], lNum);
			Light *light = scene->lights[lNum];
			// Sample ray leaving light source
			RayDifferential ray;
			float pdf;
			Spectrum alpha =
				light->Sample_L(scene, lightSamp0[2*sampOffset],
						lightSamp0[2*sampOffset+1],
						lightSamp1[2*sampOffset],
						lightSamp1[2*sampOffset+1],
						&ray, &pdf);
			if (pdf == 0.f || alpha.Black()) continue;
			alpha /= pdf * lightPdf;
//			fprintf(stderr, "initial alpha %f, light # %d\n", alpha.y(), lNum);
			Intersection isect;
			int nIntersections = 0;
			while (scene->Intersect(ray, &isect) && !alpha.Black()) {
				++nIntersections;
				alpha *= scene->Transmittance(ray);
				Vector wo = -ray.d;
				BSDF *bsdf = isect.GetBSDF(ray);
				// Create virtual light at ray intersection point
				static StatsCounter vls("IGI Integrator", "Virtual Lights Created"); //NOBOOK
				++vls; //NOBOOK
				Spectrum Le = alpha * bsdf->rho(wo) / M_PI;
//				fprintf(stderr, "\tmade light with le y %f\n", Le.y());
				virtualLights[s].push_back(VirtualLight(isect.dg.p, isect.dg.nn, Le));
				// Sample new ray direction and update weight
				Vector wi;
				float pdf;
				BxDFType flags;
				Spectrum fr = bsdf->Sample_f(wo, &wi, RandomFloat(),
								 RandomFloat(), RandomFloat(),
								 &pdf, BSDF_ALL, &flags);
				if (fr.Black() || pdf == 0.f)
					break;
				Spectrum anew = alpha * fr * AbsDot(wi, bsdf->dgShading.nn) / pdf;
				float r = anew.y() / alpha.y();
//				fprintf(stderr, "\tr = %f\n", r);
				if (RandomFloat() > r)
					break;
				alpha = anew / r;
//				fprintf(stderr, "\tnew alpha %f\n", alpha.y());
				ray = RayDifferential(isect.dg.p, wi);
			}
			BSDF::FreeAll();
		}
	}
	delete[] lightNum; // NOBOOK
	delete[] lightSamp0; // NOBOOK
	delete[] lightSamp1; // NOBOOK
}
Spectrum IGIIntegrator::Li(const Scene *scene,
		const RayDifferential &ray, const Sample *sample,
		   float *alpha) const {
	Spectrum L(0.);
	Intersection isect;
	if (scene->Intersect(ray, &isect)) {
		if (alpha) *alpha = 1.;
		Vector wo = -ray.d;
		// Compute emitted light if ray hit an area light source
		L += isect.Le(wo);
		// Evaluate BSDF at hit point
		BSDF *bsdf = isect.GetBSDF(ray);
		const Point &p = bsdf->dgShading.p;
		const Normal &n = bsdf->dgShading.nn;
		L += UniformSampleAllLights(scene, p, n,
					    wo, bsdf, sample,
					    lightSampleOffset, bsdfSampleOffset,
					    bsdfComponentOffset);
		// Compute indirect illumination with virtual lights
		u_int lSet = min(u_int(sample->oneD[vlSetOffset][0] * nLightSets),
		                 nLightSets-1);
		for (u_int i = 0; i < virtualLights[lSet].size(); ++i) {
			const VirtualLight &vl = virtualLights[lSet][i];
			// Add contribution from _VirtualLight_ _vl_
			// Ignore light if it's too close to current point
			float d2 = DistanceSquared(p, vl.p);
			//if (d2 < .8 * minDist2) continue;
			float distScale = SmoothStep(.8 * minDist2, 1.2 * minDist2, d2);
			// Compute virtual light's tentative contribution _Llight_
			Vector wi = Normalize(vl.p - p);
			Spectrum f = distScale * bsdf->f(wo, wi);
			if (f.Black()) continue;
			float G = AbsDot(wi, n) * AbsDot(wi, vl.n) / d2;
			Spectrum Llight = indirectScale * f * G * vl.Le /
				virtualLights[lSet].size();
			Llight *= scene->Transmittance(Ray(p, vl.p - p));
			// Possibly skip shadow ray with Russian roulette
			if (Llight.y() < rrThreshold) {
				float continueProbability = .1f;
				if (RandomFloat() > continueProbability)
					continue;
				Llight /= continueProbability;
			}
			static StatsCounter vlsr("IGI Integrator", "Shadow Rays to Virtual Lights"); //NOBOOK
			++vlsr; //NOBOOK
			if (!scene->IntersectP(Ray(p, vl.p - p, RAY_EPSILON,
					1.f - RAY_EPSILON)))
				L += Llight;
		}
		// Trace rays for specular reflection and refraction
		if (specularDepth++ < maxSpecularDepth) {
			Vector wi;
			// Trace rays for specular reflection and refraction
			Spectrum f = bsdf->Sample_f(wo, &wi,
						BxDFType(BSDF_REFLECTION | BSDF_SPECULAR));
			if (!f.Black()) {
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
			if (!f.Black()) {
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
		--specularDepth;
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
extern "C" DLLEXPORT SurfaceIntegrator *CreateSurfaceIntegrator(const ParamSet &params)
{
	int nLightPaths = params.FindOneInt("nlights", 64);
	int nLightSets = params.FindOneInt("nsets", 4);
	float minDist = params.FindOneFloat("mindist", .1f);
	float rrThresh = params.FindOneFloat("rrthreshold", .05f);
	float is = params.FindOneFloat("indirectscale", 1.f);
	return new IGIIntegrator(nLightPaths, nLightSets, minDist, rrThresh, is);
}
