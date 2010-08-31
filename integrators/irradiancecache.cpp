
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

// irradiancecache.cpp*
#include "pbrt.h"
#include "transport.h"
#include "scene.h"
#include "mc.h"
#include "octree.h"
// IrradianceCache Forward Declarations
struct IrradianceSample;
struct IrradProcess;
// IrradianceCache Declarations
class IrradianceCache : public SurfaceIntegrator {
public:
	// IrradianceCache Public Methods
	IrradianceCache(int maxspec, int maxind, float maxerr, int nsamples);
	~IrradianceCache();
	Spectrum Li(const Scene *scene, const RayDifferential &ray, const Sample *sample, float *alpha) const;
	void RequestSamples(Sample *sample, const Scene *scene);
	void Preprocess(const Scene *);
private:
	// IrradianceCache Data
	float maxError;
	int nSamples;
	int maxSpecularDepth, maxIndirectDepth;
	mutable int specularDepth;
	// Declare sample parameters for light source sampling
	int *lightSampleOffset, lightNumOffset;
	int *bsdfSampleOffset, *bsdfComponentOffset;
	mutable Octree<IrradianceSample, IrradProcess> *octree;
	// IrradianceCache Private Methods
	Spectrum IndirectLo(const Point &p, const Normal &n,
		const Vector &wo, BSDF *bsdf, BxDFType flags,
		const Sample *sample, const Scene *scene) const;
	bool InterpolateIrradiance(const Scene *scene,
			const Point &p, const Normal &n, Spectrum *E) const;
};
struct IrradianceSample {
	// IrradianceSample Constructor
	IrradianceSample() { }
	IrradianceSample(const Spectrum &e, const Point &P, const Normal &N,
		float md) : E(e), n(N), p(P) {
		maxDist = md;
	}
	Spectrum E;
	Normal n;
	Point p;
	float maxDist;
};
struct IrradProcess {
	// IrradProcess Public Methods
	IrradProcess(const Normal &N, float me) {
		n = N;
		maxError = me;
		nFound = samplesChecked = 0;
		sumWt = 0.;
		E = 0.;
	}
	void operator()(const Point &P, const IrradianceSample &sample) const;
	bool Successful() {
		return (sumWt > 0. && nFound > 0);
	}
	Spectrum GetIrradiance() const { return E / sumWt; }
	Normal n;
	float maxError;
	mutable int nFound, samplesChecked;
	mutable float sumWt;
	mutable Spectrum E;
};
// IrradianceCache Method Definitions
IrradianceCache::IrradianceCache(int maxspec, int maxind,
		float maxerr, int ns) {
	maxError = maxerr;
	nSamples = ns;
	maxSpecularDepth = maxspec;
	maxIndirectDepth = maxind;
	specularDepth = 0;
}
void IrradianceCache::RequestSamples(Sample *sample,
		const Scene *scene) {
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
Spectrum IrradianceCache::Li(const Scene *scene, const RayDifferential &ray,
		const Sample *sample, float *alpha) const {
	Intersection isect;
	Spectrum L(0.);
	if (scene->Intersect(ray, &isect)) {
		if (alpha) *alpha = 1.;
		// Evaluate BSDF at hit point
		BSDF *bsdf = isect.GetBSDF(ray);
		Vector wo = -ray.d;
		const Point &p = bsdf->dgShading.p;
		const Normal &n = bsdf->dgShading.nn;
		// Compute direct lighting for irradiance cache
		L += isect.Le(wo);
		L += UniformSampleAllLights(scene, p, n, wo, bsdf, sample,
			lightSampleOffset, bsdfSampleOffset,
			bsdfComponentOffset);
		// Compute indirect lighting for irradiance cache
		if (specularDepth++ < maxSpecularDepth) {
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
		--specularDepth;
		// Estimate indirect lighting with irradiance cache
		Normal ng = isect.dg.nn;
		if (Dot(wo, ng) < 0.f) ng = -ng;
		BxDFType flags = BxDFType(BSDF_REFLECTION |
		                          BSDF_DIFFUSE |
								  BSDF_GLOSSY);
		L += IndirectLo(p, ng, wo, bsdf, flags, sample, scene);
		flags = BxDFType(BSDF_TRANSMISSION |
		                 BSDF_DIFFUSE |
						 BSDF_GLOSSY);
		L += IndirectLo(p, -ng, wo, bsdf, flags, sample, scene);
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
Spectrum IrradianceCache::IndirectLo(const Point &p,
		const Normal &n, const Vector &wo, BSDF *bsdf,
		BxDFType flags, const Sample *sample,
		const Scene *scene) const {
	if (bsdf->NumComponents(flags) == 0)
		return Spectrum(0.);
	Spectrum E;
	if (!InterpolateIrradiance(scene, p, n, &E)) {
		// Compute irradiance at current point
		u_int scramble[2] = { RandomUInt(), RandomUInt() };
		float sumInvDists = 0.;
		for (int i = 0; i < nSamples; ++i) {
			// Trace ray to sample radiance for irradiance estimate
			// Update irradiance statistics for rays traced
			static StatsCounter nIrradiancePaths("Irradiance Cache",
				"Paths followed for irradiance estimates");
			++nIrradiancePaths;
			float u[2];
			Sample02(i, scramble, u);
			Vector w = CosineSampleHemisphere(u[0], u[1]);
			RayDifferential r(p, bsdf->LocalToWorld(w));
			if (Dot(r.d, n) < 0) r.d = -r.d;
			Spectrum L(0.);
			// Do path tracing to compute radiance along ray for estimate
			{
			// Declare common path integration variables
			Spectrum pathThroughput = 1.;
			RayDifferential ray(r);
			bool specularBounce = false;
			for (int pathLength = 0; ; ++pathLength) {
				// Find next vertex of path
				Intersection isect;
				if (!scene->Intersect(ray, &isect))
					break;
				if (pathLength == 0)
					r.maxt = ray.maxt;
				pathThroughput *= scene->Transmittance(ray);
				// Possibly add emitted light at path vertex
				if (specularBounce)
					L += pathThroughput * isect.Le(-ray.d);
				// Evaluate BSDF at hit point
				BSDF *bsdf = isect.GetBSDF(ray);
				// Sample illumination from lights to find path contribution
				const Point &p = bsdf->dgShading.p;
				const Normal &n = bsdf->dgShading.nn;
				Vector wo = -ray.d;
				L += pathThroughput *
					UniformSampleOneLight(scene, p, n, wo, bsdf, sample);
				if (pathLength+1 == maxIndirectDepth) break;
				// Sample BSDF to get new path direction
				// Get random numbers for sampling new direction, \mono{bs1}, \mono{bs2}, and \mono{bcs}
				float bs1 = RandomFloat(), bs2 = RandomFloat(), bcs = RandomFloat();
				Vector wi;
				float pdf;
				BxDFType flags;
				Spectrum f = bsdf->Sample_f(wo, &wi, bs1, bs2, bcs,
					&pdf, BSDF_ALL, &flags);
				if (f.Black() || pdf == 0.)
					break;
				specularBounce = (flags & BSDF_SPECULAR) != 0;
				pathThroughput *= f * AbsDot(wi, n) / pdf;
				ray = RayDifferential(p, wi);
				// Possibly terminate the path
				if (pathLength > 3) {
					float continueProbability = .5f;
					if (RandomFloat() > continueProbability)
						break;
					pathThroughput /= continueProbability;
				}
			}
			}
			E += L;
			float dist = r.maxt * r.d.Length();
			sumInvDists += 1.f / dist;
		}
		E *= M_PI / float(nSamples);
		// Add computed irradiance value to cache
		// Update statistics for new irradiance sample
		static StatsCounter nSamplesComputed("Irradiance Cache",
			"Irradiance estimates computed");
		++nSamplesComputed;
		// Compute bounding box of irradiance sample's contribution region
		static float minMaxDist =
			.001f * powf(scene->WorldBound().Volume(), 1.f/3.f);
		static float maxMaxDist =
			.125f * powf(scene->WorldBound().Volume(), 1.f/3.f);
		float maxDist = nSamples / sumInvDists;
		if (minMaxDist > 0.f)
			maxDist = Clamp(maxDist, minMaxDist, maxMaxDist);
		maxDist *= maxError;
		BBox sampleExtent(p);
		sampleExtent.Expand(maxDist);
		octree->Add(IrradianceSample(E, p, n, maxDist),
			sampleExtent);
	}
	return INV_PI * bsdf->rho(wo, flags) * E;
}
void IrradianceCache::Preprocess(const Scene *scene) {
	BBox wb = scene->WorldBound();
	Vector delta = .01f * (wb.pMax - wb.pMin);
	wb.pMin -= delta;
	wb.pMax += delta;
	octree = new Octree<IrradianceSample, IrradProcess>(wb);
}
IrradianceCache::~IrradianceCache() {
	delete octree;
}
bool
IrradianceCache::InterpolateIrradiance(const Scene *scene,
		const Point &p, const Normal &n, Spectrum *E) const {
	if (!octree) return false;
	IrradProcess proc(n, maxError);
	octree->Lookup(p, proc);
	// Update irradiance cache lookup statistics
	static StatsPercentage nSuccessfulLookups("Irradiance Cache",
		"Successful irradiance cache lookups");
	static StatsRatio nSamplesFound("Irradiance Cache",
		"Irradiance samples found per successful lookup");
	static StatsRatio nSamplesChecked("Irradiance Cache",
		"Irradiance samples checked per lookup");
	nSuccessfulLookups.Add(proc.Successful() ? 1 : 0, 1);
	nSamplesFound.Add(proc.nFound, 1);
	nSamplesChecked.Add(proc.samplesChecked, 1);
	if (!proc.Successful()) return false;
	*E = proc.GetIrradiance();
	return true;
}
void IrradProcess::operator()(const Point &p,
		const IrradianceSample &sample) const {
	++samplesChecked;
	// Skip irradiance sample if surface normals are too different
	if (Dot(n, sample.n) < 0.01f)
		return;
	// Skip irradiance sample if it's too far from the sample point
	float d2 = DistanceSquared(p, sample.p);
	if (d2 > sample.maxDist * sample.maxDist)
		return;
	// Skip irradiance sample if it's in front of point being shaded
	Normal navg = sample.n + n;
	if (Dot(p - sample.p, navg) < -.01f)
		return;
	// Compute estimate error term and possibly use sample
	float err = sqrtf(d2) / (sample.maxDist * Dot(n, sample.n));
	if (err < 1.) {
		++nFound;
		float wt = (1.f - err) * (1.f - err);
		E += wt * sample.E;
		sumWt += wt;
	}
}
extern "C" DLLEXPORT SurfaceIntegrator *CreateSurfaceIntegrator(const ParamSet &params) {
	float maxError = params.FindOneFloat("maxerror", .2f);
	int maxSpecularDepth = params.FindOneInt("maxspeculardepth", 5);
	int maxIndirectDepth = params.FindOneInt("maxindirectdepth", 3);
	int nSamples = params.FindOneInt("nsamples", 4096);
	return new IrradianceCache(maxSpecularDepth, maxIndirectDepth,
		maxError, nSamples);
}
