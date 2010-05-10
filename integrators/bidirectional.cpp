
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

// bidirectional.cpp*
#include "pbrt.h"
#include "transport.h"
#include "scene.h"
#include "mc.h"
// Bidirectional Local Declarations
struct BidirVertex;
class BidirIntegrator : public SurfaceIntegrator {
public:
	// BidirIntegrator Public Methods
	Spectrum Li(const Scene *scene, const RayDifferential &ray, const Sample *sample, float *alpha) const;
	void RequestSamples(Sample *sample, const Scene *scene);
private:
	// BidirIntegrator Private Methods
	int generatePath(const Scene *scene, const Ray &r, const Sample *sample,
		const int *bsdfOffset, const int *bsdfCompOffset,
		BidirVertex *vertices, int maxVerts) const;
	float weightPath(BidirVertex *eye, int nEye, BidirVertex *light, int nLight) const;
	Spectrum evalPath(const Scene *scene, BidirVertex *eye, int nEye,
	BidirVertex *light, int nLight) const;
	static float G(const BidirVertex &v0, const BidirVertex &v1);
	static bool visible(const Scene *scene, const Point &P0, const Point &P1);
	// BidirIntegrator Data
	#define MAX_VERTS 4
	int eyeBSDFOffset[MAX_VERTS], eyeBSDFCompOffset[MAX_VERTS];
	int lightBSDFOffset[MAX_VERTS], lightBSDFCompOffset[MAX_VERTS];
	int directLightOffset[MAX_VERTS], directLightNumOffset[MAX_VERTS];
	int directBSDFOffset[MAX_VERTS], directBSDFCompOffset[MAX_VERTS];
	int lightNumOffset, lightPosOffset, lightDirOffset;
};
struct BidirVertex {
	BidirVertex() { bsdfWeight = dAWeight = 0.; rrWeight = 1.;
		flags = BxDFType(0); bsdf = NULL; }
	BSDF *bsdf;
	Point p;
	Normal ng, ns;
	Vector wi, wo;
	float bsdfWeight, dAWeight, rrWeight;
	BxDFType flags;
};
// Bidirectional Method Definitions
void BidirIntegrator::RequestSamples(Sample *sample, const Scene *scene) {
	for (int i = 0; i < MAX_VERTS; ++i) {
		eyeBSDFOffset[i] = sample->Add2D(1);
		eyeBSDFCompOffset[i] = sample->Add1D(1);
		lightBSDFOffset[i] = sample->Add2D(1);
		lightBSDFCompOffset[i] = sample->Add1D(1);
		directLightOffset[i] = sample->Add2D(1);
		directLightNumOffset[i] = sample->Add1D(1);
		directBSDFOffset[i] = sample->Add2D(1);
		directBSDFCompOffset[i]  = sample->Add1D(1);
	}
	lightNumOffset = sample->Add1D(1);
	lightPosOffset = sample->Add2D(1);
	lightDirOffset = sample->Add2D(1);
}
Spectrum BidirIntegrator::Li(const Scene *scene,
		const RayDifferential &ray,
		const Sample *sample, float *alpha) const {
	Spectrum L(0.);
	// Generate eye and light sub-paths
	BidirVertex eyePath[MAX_VERTS], lightPath[MAX_VERTS];
	int nEye = generatePath(scene, ray, sample, eyeBSDFOffset,
		eyeBSDFCompOffset, eyePath, MAX_VERTS);
	if (nEye == 0) {
		*alpha = 0.;
		return L;
	}
	*alpha = 1;
	// Choose light for bidirectional path
	int lightNum = Floor2Int(sample->oneD[lightNumOffset][0] *
		scene->lights.size());
	lightNum = min(lightNum, (int)scene->lights.size() - 1);
	Light *light = scene->lights[lightNum];
	float lightWeight = float(scene->lights.size());
	// Sample ray from light source to start light path
	Ray lightRay;
	float lightPdf;
	float u[4];
	u[0] = sample->twoD[lightPosOffset][0];
	u[1] = sample->twoD[lightPosOffset][1];
	u[2] = sample->twoD[lightDirOffset][0];
	u[3] = sample->twoD[lightDirOffset][1];
	Spectrum Le = light->Sample_L(scene, u[0], u[1], u[2], u[3],
		&lightRay, &lightPdf);
	if (lightPdf == 0.) return 0.f;
	Le = lightWeight / lightPdf;
	int nLight = generatePath(scene, lightRay, sample, lightBSDFOffset,
		lightBSDFCompOffset, lightPath, MAX_VERTS);
	// Connect bidirectional path prefixes and evaluate throughput
	Spectrum directWt(1.0);
	for (int i = 1; i <= nEye; ++i) {
		// Handle direct lighting for bidirectional integrator
		directWt /= eyePath[i-1].rrWeight;
		L += directWt *
			UniformSampleOneLight(scene, eyePath[i-1].p, eyePath[i-1].ng, eyePath[i-1].wi,
			eyePath[i-1].bsdf, sample, directLightOffset[i-1], directLightNumOffset[i-1],
			directBSDFOffset[i-1], directBSDFCompOffset[i-1]) /
			weightPath(eyePath, i, lightPath, 0);
		directWt *= eyePath[i-1].bsdf->f(eyePath[i-1].wi, eyePath[i-1].wo) *
			AbsDot(eyePath[i-1].wo, eyePath[i-1].ng) /
			eyePath[i-1].bsdfWeight;
		for (int j = 1; j <= nLight; ++j)
			L += Le * evalPath(scene, eyePath, i, lightPath, j) /
				weightPath(eyePath, i, lightPath, j);
	}
	return L;
}
int BidirIntegrator::generatePath(const Scene *scene, const Ray &r,
		const Sample *sample, const int *bsdfOffset,
		const int *bsdfCompOffset,
		BidirVertex *vertices, int maxVerts) const {
	int nVerts = 0;
	RayDifferential ray(r.o, r.d);
	while (nVerts < maxVerts) {
		// Find next vertex in path and initialize _vertices_
		Intersection isect;
		if (!scene->Intersect(ray, &isect))
			break;
		BidirVertex &v = vertices[nVerts];
		v.bsdf = isect.GetBSDF(ray); // do before Ns is set!
		v.p = isect.dg.p;
		v.ng = isect.dg.nn;
		v.ns = v.bsdf->dgShading.nn;
		v.wi = -ray.d;
		++nVerts;
		// Possibly terminate bidirectional path sampling
		if (nVerts > 2) {
			float rrProb = .2f;
			if (RandomFloat() > rrProb)
				break;
			v.rrWeight = 1.f / rrProb;
		}
		// Initialize _ray_ for next segment of path
		float u1 = sample->twoD[bsdfOffset[nVerts-1]][0];
		float u2 = sample->twoD[bsdfOffset[nVerts-1]][1];
		float u3 = sample->oneD[bsdfCompOffset[nVerts-1]][0];
		Spectrum fr = v.bsdf->Sample_f(v.wi, &v.wo, u1, u2, u3,
			 &v.bsdfWeight, BSDF_ALL, &v.flags);
		if (fr.Black() && v.bsdfWeight == 0.f)
			break;
		ray = RayDifferential(v.p, v.wo);
	}
	// Initialize additional values in _vertices_
	for (int i = 0; i < nVerts-1; ++i)
		vertices[i].dAWeight = vertices[i].bsdfWeight *
			AbsDot(-vertices[i].wo, vertices[i+1].ng) /
			DistanceSquared(vertices[i].p, vertices[i+1].p);
	return nVerts;
}

float BidirIntegrator::weightPath(BidirVertex *eye, int nEye,
		BidirVertex *light, int nLight) const {
	return float(nEye + nLight);
}
Spectrum BidirIntegrator::evalPath(const Scene *scene, BidirVertex *eye, int nEye,
		BidirVertex *light, int nLight) const {
	Spectrum L(1.);
	for (int i = 0; i < nEye-1; ++i)
		L *= eye[i].bsdf->f(eye[i].wi, eye[i].wo) *
			AbsDot(eye[i].wo, eye[i].ng) /
			(eye[i].bsdfWeight * eye[i].rrWeight);
	Vector w = light[nLight-1].p - eye[nEye-1].p;
	L *= eye[nEye-1].bsdf->f(eye[nEye-1].wi, w) *
		G(eye[nEye-1], light[nLight-1]) *
		light[nLight-1].bsdf->f(-w, light[nLight-1].wi) /
		(eye[nEye-1].rrWeight * light[nLight-1].rrWeight);
	for (int i = nLight-2; i >= 0; --i)
		L *= light[i].bsdf->f(light[i].wi, light[i].wo) *
			AbsDot(light[i].wo, light[i].ng) /
			(light[i].bsdfWeight * light[i].rrWeight);
	if (L.Black())
		return L;
	if (!visible(scene, eye[nEye-1].p, light[nLight-1].p))
		return 0.;
	return L;
}
float BidirIntegrator::G(const BidirVertex &v0, const BidirVertex &v1) {
	Vector w = Normalize(v1.p - v0.p);
	return AbsDot(v0.ng, w) * AbsDot(v1.ng, -w) /
		DistanceSquared(v0.p, v1.p);
}
bool BidirIntegrator::visible(const Scene *scene, const Point &P0,
		const Point &P1) {
	Ray ray(P0, P1-P0, RAY_EPSILON, 1.f - RAY_EPSILON);
	return !scene->IntersectP(ray);
}
extern "C" DLLEXPORT SurfaceIntegrator *CreateSurfaceIntegrator(const ParamSet &params) {
	return new BidirIntegrator;
}
