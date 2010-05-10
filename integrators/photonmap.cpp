
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

// photonmap.cpp*
#include "pbrt.h"
#include "transport.h"
#include "scene.h"
#include "mc.h"
#include "kdtree.h"
#include "sampling.h"
// Photonmap Local Declarations
struct Photon;
struct ClosePhoton;
struct PhotonProcess;
class PhotonIntegrator : public SurfaceIntegrator {
public:
	// PhotonIntegrator Public Methods
	PhotonIntegrator(int ncaus, int ndir, int nindir, int nLookup, int mdepth,
		float maxdist, bool finalGather, int gatherSamples,
		bool directWithPhotons);
	~PhotonIntegrator();
	Spectrum Li(const Scene *scene, const RayDifferential &ray,
		const Sample *sample, float *alpha) const;
	void RequestSamples(Sample *sample, const Scene *scene);
	void Preprocess(const Scene *);
private:
	// PhotonIntegrator Private Methods
	static inline bool unsuccessful(int needed, int found, int shot) {
		return (found < needed &&
		        (found == 0 || found < shot / 1024));
	}
	static Spectrum LPhoton(KdTree<Photon, PhotonProcess> *map,
		int nPaths, int nLookup, BSDF *bsdf, const Intersection &isect,
		const Vector &w, float maxDistSquared);
	// PhotonIntegrator Private Data
	u_int nCausticPhotons, nIndirectPhotons, nDirectPhotons;
	u_int nLookup;
	mutable int specularDepth;
	int maxSpecularDepth;
	float maxDistSquared;
	bool directWithPhotons, finalGather;
	int gatherSamples;
	// Declare sample parameters for light source sampling
	int *lightSampleOffset, lightNumOffset;
	int *bsdfSampleOffset, *bsdfComponentOffset;
	int gatherSampleOffset, gatherComponentOffset;
	int nCausticPaths, nDirectPaths, nIndirectPaths;
	mutable KdTree<Photon, PhotonProcess> *causticMap;
	mutable KdTree<Photon, PhotonProcess> *directMap;
	mutable KdTree<Photon, PhotonProcess> *indirectMap;
};
struct Photon {
	// Photon Constructor
	Photon(const Point &pp, const Spectrum &wt, const Vector &w)
		: p(pp), alpha(wt), wi(w) {
	}
	Photon() { }
	Point p;
	Spectrum alpha;
	Vector wi;
};
struct PhotonProcess {
	// PhotonProcess Public Methods
	PhotonProcess(u_int mp, const Point &p);
	void operator()(const Photon &photon, float dist2, float &maxDistSquared) const;
	const Point &p;
	ClosePhoton *photons;
	u_int nLookup;
	mutable u_int foundPhotons;
};
struct ClosePhoton {
	ClosePhoton(const Photon *p = NULL,
	            float md2 = INFINITY) {
		photon = p;
		distanceSquared = md2;
	}
	bool operator<(const ClosePhoton &p2) const {
		return distanceSquared == p2.distanceSquared ? (photon < p2.photon) :
			distanceSquared < p2.distanceSquared;
	}
	const Photon *photon;
	float distanceSquared;
};
// Photonmap Method Definitions
PhotonIntegrator::PhotonIntegrator(int ncaus, int ndir, int nind,
		int nl,	int mdepth, float mdist, bool fg,
		int gs, bool dp) {
	nCausticPhotons = ncaus;
	nIndirectPhotons = nind;
	nDirectPhotons = ndir;
	nLookup = nl;
	maxDistSquared = mdist * mdist;
	maxSpecularDepth = mdepth;
	causticMap = directMap = indirectMap = NULL;
	specularDepth = 0;
	finalGather = fg;
	gatherSamples = gs;
	directWithPhotons = dp;
}
PhotonIntegrator::~PhotonIntegrator() {
	delete causticMap;
	delete directMap;
	delete indirectMap;
}
void PhotonIntegrator::RequestSamples(Sample *sample,
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
	if (finalGather) {
		gatherSamples = scene->sampler->RoundSize(gatherSamples);
		gatherSampleOffset = sample->Add2D(gatherSamples);
		gatherComponentOffset = sample->Add1D(gatherSamples);
	}
}
void PhotonIntegrator::Preprocess(const Scene *scene) {
	if (scene->lights.size() == 0) return;
	ProgressReporter progress(nCausticPhotons+nDirectPhotons+ // NOBOOK
		nIndirectPhotons, "Shooting photons"); // NOBOOK
	vector<Photon> causticPhotons;
	vector<Photon> directPhotons;
	vector<Photon> indirectPhotons;
	causticPhotons.reserve(nCausticPhotons); // NOBOOK
	directPhotons.reserve(nDirectPhotons); // NOBOOK
	indirectPhotons.reserve(nIndirectPhotons); // NOBOOK
	// Initialize photon shooting statistics
	static StatsCounter nshot("Photon Map",
		"Number of photons shot from lights");
	bool causticDone = (nCausticPhotons == 0);
	bool directDone = (nDirectPhotons == 0);
	bool indirectDone = (nIndirectPhotons == 0);
	while (!causticDone || !directDone || !indirectDone) {
		++nshot;
		// Give up if we're not storing enough photons
		if (nshot > 500000 &&
			(unsuccessful(nCausticPhotons,
			              causticPhotons.size(),
						  nshot) ||
			 unsuccessful(nDirectPhotons,
			              directPhotons.size(),
						  nshot) ||
			 unsuccessful(nIndirectPhotons,
			              indirectPhotons.size(),
						  nshot))) {
			Error("Unable to store enough photons.  Giving up.\n");
			return;
		}
		// Trace a photon path and store contribution
		// Choose 4D sample values for photon
		float u[4];
		u[0] = (float)RadicalInverse((int)nshot+1, 2);
		u[1] = (float)RadicalInverse((int)nshot+1, 3);
		u[2] = (float)RadicalInverse((int)nshot+1, 5);
		u[3] = (float)RadicalInverse((int)nshot+1, 7);
		// Choose light to shoot photon from
		int nLights = int(scene->lights.size());
		int lightNum =
			min(Floor2Int(nLights * (float)RadicalInverse((int)nshot+1, 11)),
			nLights-1);
		Light *light = scene->lights[lightNum];
		float lightPdf = 1.f / nLights;
		// Generate _photonRay_ from light source and initialize _alpha_
		RayDifferential photonRay;
		float pdf;
		Spectrum alpha =
			light->Sample_L(scene, u[0], u[1], u[2], u[3],
				&photonRay, &pdf);
		if (pdf == 0.f || alpha.Black()) continue;
		alpha /= pdf * lightPdf;
		if (!alpha.Black()) {
			// Follow photon path through scene and record intersections
			bool specularPath = false;
			Intersection photonIsect;
			int nIntersections = 0;
			while (scene->Intersect(photonRay, &photonIsect)) {
				++nIntersections;
				// Handle photon/surface intersection
				alpha *= scene->Transmittance(photonRay);
				Vector wo = -photonRay.d;
				BSDF *photonBSDF = photonIsect.GetBSDF(photonRay);
				BxDFType specularType = BxDFType(BSDF_REFLECTION |
					BSDF_TRANSMISSION | BSDF_SPECULAR);
				bool hasNonSpecular = (photonBSDF->NumComponents() >
					photonBSDF->NumComponents(specularType));
				if (hasNonSpecular) {
					// Deposit photon at surface
					Photon photon(photonIsect.dg.p, alpha, wo);
					if (nIntersections == 1) {
						// Process direct lighting photon intersection
						if (!directDone) {
							directPhotons.push_back(photon);
							if (directPhotons.size() == nDirectPhotons) {
								directDone = true;
								nDirectPaths = (int)nshot;
								directMap =
									new KdTree<Photon,
											   PhotonProcess>(directPhotons);
							}
							progress.Update(); // NOBOOK
						}
					}
					else if (specularPath) {
						// Process caustic photon intersection
						if (!causticDone) {
							causticPhotons.push_back(photon);
							if (causticPhotons.size() == nCausticPhotons) {
								causticDone = true;
								nCausticPaths = (int)nshot;
								causticMap =
									new KdTree<Photon,
										       PhotonProcess>(causticPhotons);
							}
							progress.Update();
						}
					}
					else {
						// Process indirect lighting photon intersection
						if (!indirectDone) {
							indirectPhotons.push_back(photon);
							if (indirectPhotons.size() == nIndirectPhotons) {
								indirectDone = true;
								nIndirectPaths = (int)nshot;
								indirectMap =
									new KdTree<Photon,
											   PhotonProcess>(indirectPhotons);
							}
							progress.Update();
						}
					}
				}
				// Sample new photon ray direction
				Vector wi;
				float pdf;
				BxDFType flags;
				// Get random numbers for sampling outgoing photon direction
				float u1, u2, u3;
				if (nIntersections == 1) {
					u1 = (float)RadicalInverse((int)nshot+1, 13);
					u2 = (float)RadicalInverse((int)nshot+1, 17);
					u3 = (float)RadicalInverse((int)nshot+1, 19);
				}
				else {
					u1 = RandomFloat();
					u2 = RandomFloat();
					u3 = RandomFloat();
				}
				Spectrum fr = photonBSDF->Sample_f(wo, &wi, u1, u2, u3,
					&pdf, BSDF_ALL, &flags);
				if (fr.Black() || pdf == 0.f)
					break;
				specularPath = (nIntersections == 1 || specularPath) &&
					((flags & BSDF_SPECULAR) != 0);
				alpha *= fr * AbsDot(wi, photonBSDF->dgShading.nn) / pdf;
				photonRay = RayDifferential(photonIsect.dg.p, wi);
				// Possibly terminate photon path
				if (nIntersections > 3) {
					float continueProbability = .5f;
					if (RandomFloat() > continueProbability)
						break;
					alpha /= continueProbability;
				}
			}
		}
		BSDF::FreeAll();
	}
	progress.Done(); // NOBOOK
}
Spectrum PhotonIntegrator::Li(const Scene *scene,
		const RayDifferential &ray, const Sample *sample,
		float *alpha) const {
	// Compute reflected radiance with photon map
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
		// Compute direct lighting for photon map integrator
		if (directWithPhotons)
			L += LPhoton(directMap, nDirectPaths, nLookup,
				bsdf, isect, wo, maxDistSquared);
		else
			L += UniformSampleAllLights(scene, p, n,
				wo, bsdf, sample,
				lightSampleOffset, bsdfSampleOffset,
				bsdfComponentOffset);
		
		// Compute indirect lighting for photon map integrator
		L += LPhoton(causticMap, nCausticPaths, nLookup, bsdf,
			isect, wo, maxDistSquared);
		if (finalGather) {
			// Do one-bounce final gather for photon map
			Spectrum Li(0.);
			for (int i = 0; i < gatherSamples; ++i) {
				// Sample random direction for final gather ray
				Vector wi;
				float u1 = sample->twoD[gatherSampleOffset][2*i];
				float u2 = sample->twoD[gatherSampleOffset][2*i+1];
				float u3 = sample->oneD[gatherComponentOffset][i];
				float pdf;
				Spectrum fr = bsdf->Sample_f(wo, &wi, u1, u2, u3,
					&pdf, BxDFType(BSDF_ALL & (~BSDF_SPECULAR)));
				if (fr.Black() || pdf == 0.f) continue;
				RayDifferential bounceRay(p, wi);
				static StatsCounter gatherRays("Photon Map", // NOBOOK
					"Final gather rays traced"); // NOBOOK
				++gatherRays; // NOBOOK
				Intersection gatherIsect;
				if (scene->Intersect(bounceRay, &gatherIsect)) {
					// Compute exitant radiance at final gather intersection
					BSDF *gatherBSDF = gatherIsect.GetBSDF(bounceRay);
					Vector bounceWo = -bounceRay.d;
					Spectrum Lindir =
						LPhoton(directMap, nDirectPaths, nLookup,
							gatherBSDF, gatherIsect, bounceWo, maxDistSquared) +
						LPhoton(indirectMap, nIndirectPaths, nLookup,
							gatherBSDF, gatherIsect, bounceWo, maxDistSquared) +
						LPhoton(causticMap, nCausticPaths, nLookup,
							gatherBSDF, gatherIsect, bounceWo, maxDistSquared);
					Lindir *= scene->Transmittance(bounceRay);
					Li += fr * Lindir * AbsDot(wi, n) / pdf;
				}
			}
			L += Li / float(gatherSamples);
		}
		else
			L += LPhoton(indirectMap, nIndirectPaths, nLookup,
				bsdf, isect, wo, maxDistSquared);
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
Spectrum PhotonIntegrator::LPhoton(
		KdTree<Photon, PhotonProcess> *map,
		int nPaths, int nLookup, BSDF *bsdf,
		const Intersection &isect, const Vector &wo,
		float maxDistSquared) {
	Spectrum L(0.);
	if (!map) return L;
	BxDFType nonSpecular = BxDFType(BSDF_REFLECTION |
		BSDF_TRANSMISSION | BSDF_DIFFUSE | BSDF_GLOSSY);
	if (bsdf->NumComponents(nonSpecular) == 0)
		return L;
	static StatsCounter lookups("Photon Map", "Total lookups"); // NOBOOK
	// Initialize _PhotonProcess_ object, _proc_, for photon map lookups
	PhotonProcess proc(nLookup, isect.dg.p);
	proc.photons =
		(ClosePhoton *)alloca(nLookup * sizeof(ClosePhoton));
	// Do photon map lookup
	++lookups;  // NOBOOK
	map->Lookup(isect.dg.p, proc, maxDistSquared);
	// Accumulate light from nearby photons
	static StatsRatio foundRate("Photon Map", "Photons found per lookup"); // NOBOOK
	foundRate.Add(proc.foundPhotons, 1); // NOBOOK
	float scale = 1.f / (float(nPaths) * maxDistSquared * M_PI);
	// Estimate reflected light from photons
	ClosePhoton *photons = proc.photons;
	int nFound = proc.foundPhotons;
	Normal Nf = Dot(wo, bsdf->dgShading.nn) < 0 ? -bsdf->dgShading.nn :
		bsdf->dgShading.nn;
	if (bsdf->NumComponents(BxDFType(BSDF_REFLECTION |
			BSDF_TRANSMISSION | BSDF_GLOSSY)) > 0) {
		// Compute exitant radiance from photons for glossy surface
		for (int i = 0; i < nFound; ++i) {
			BxDFType flag = Dot(Nf, photons[i].photon->wi) > 0.f ?
				BSDF_ALL_REFLECTION : BSDF_ALL_TRANSMISSION;
			L += bsdf->f(wo, photons[i].photon->wi,	flag) *
				(scale * photons[i].photon->alpha);
		}
	}
	else {
		// Compute exitant radiance from photons for diffuse surface
		Spectrum Lr(0.), Lt(0.);
		for (int i = 0; i < nFound; ++i)
			if (Dot(Nf, photons[i].photon->wi) > 0.f)
				Lr += photons[i].photon->alpha;
			else
				Lt += photons[i].photon->alpha;
		L += (scale * INV_PI) * (Lr * bsdf->rho(wo, BSDF_ALL_REFLECTION) +
			Lt * bsdf->rho(wo, BSDF_ALL_TRANSMISSION));
	}
	return L;
}
PhotonProcess::PhotonProcess(u_int mp, const Point &P)
	: p(P) {
	photons = 0;
	nLookup = mp;
	foundPhotons = 0;
}
void PhotonProcess::operator()(const Photon &photon,
		float distSquared, float &maxDistSquared) const {
	static StatsPercentage discarded("Photon Map", "Discarded photons"); // NOBOOK
	discarded.Add(0, 1); // NOBOOK
	if (foundPhotons < nLookup) {
		// Add photon to unordered array of photons
		photons[foundPhotons++] = ClosePhoton(&photon, distSquared);
		if (foundPhotons == nLookup) {
			std::make_heap(&photons[0], &photons[nLookup]);
			maxDistSquared = photons[0].distanceSquared;
		}
	}
	else {
		// Remove most distant photon from heap and add new photon
		discarded.Add(1, 0); // NOBOOK
		std::pop_heap(&photons[0], &photons[nLookup]);
		photons[nLookup-1] = ClosePhoton(&photon, distSquared);
		std::push_heap(&photons[0], &photons[nLookup]);
		maxDistSquared = photons[0].distanceSquared;
	}
}
extern "C" DLLEXPORT SurfaceIntegrator *CreateSurfaceIntegrator(const ParamSet &params) {
	int nCaustic = params.FindOneInt("causticphotons", 20000);
	int nDirect = params.FindOneInt("directphotons", 100000);
	int nIndirect = params.FindOneInt("indirectphotons", 100000);
	int nUsed = params.FindOneInt("nused", 50);
	int maxDepth = params.FindOneInt("maxdepth", 5);
	bool finalGather = params.FindOneBool("finalgather", true);
	bool directPhotons = params.FindOneBool("directwithphotons", false);
	int gatherSamples = params.FindOneInt("finalgathersamples", 32);
	float maxDist = params.FindOneFloat("maxdist", .1f);
	return new PhotonIntegrator(nCaustic, nDirect, nIndirect,
		nUsed, maxDepth, maxDist, finalGather, gatherSamples,
		directPhotons);
}
