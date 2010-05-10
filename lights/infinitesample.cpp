
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

// infinitesample.cpp*
#include "pbrt.h"
#include "light.h"
#include "texture.h"
#include "shape.h"
#include "scene.h"
#include "mipmap.h"
// Utility Classes and Functions
struct Distribution1D {
	// Distribution1D Methods
	Distribution1D(float *f, int n) {
		func = new float[n];
		cdf = new float[n+1];
		count = n;
		memcpy(func, f, n*sizeof(float));
		ComputeStep1dCDF(func, n, &funcInt, cdf);
		invFuncInt = 1.f / funcInt;
		invCount = 1.f / count;
	}
	float Sample(float u, float *pdf) {
		// Find surrounding cdf segments
		float *ptr = std::upper_bound(cdf, cdf+count+1, u);
		int offset = max(0, (int) (ptr-cdf-1));
		Assert(u >= cdf[offset] && u < cdf[offset+1]);
		// Return offset along current cdf segment
		u = (u - cdf[offset]) / (cdf[offset+1] - cdf[offset]);
		*pdf = func[offset] * invFuncInt;
		return offset + u;
	}
	// Distribution1D Data
	float *func, *cdf;
	float funcInt, invFuncInt, invCount;
	int count;
};
// InfiniteAreaLightIS Definitions
class InfiniteAreaLightIS : public Light {
public:
	// InfiniteAreaLightIS Public Methods
	InfiniteAreaLightIS(const Transform &light2world,	const Spectrum &power, int ns,
			  const string &texmap);
	~InfiniteAreaLightIS();
	Spectrum Power(const Scene *scene) const {
		Point worldCenter;
		float worldRadius;
		scene->WorldBound().BoundingSphere(&worldCenter,
			&worldRadius);
		return Lbase * radianceMap->Lookup(.5f, .5f, .5f) *
			M_PI * worldRadius * worldRadius;
	}
	bool IsDeltaLight() const { return false; }
	Spectrum Le(const RayDifferential &r) const;
	Spectrum Sample_L(const Point &p, float u1, float u2, Vector *wi, float *pdf,
		VisibilityTester *visibility) const;
	Spectrum Sample_L(const Scene *scene, float u1, float u2,
			float u3, float u4, Ray *ray, float *pdf) const;
	float Pdf(const Point &, const Vector &) const;
	Spectrum Sample_L(const Point &P, Vector *w, VisibilityTester *visibility) const;

private:
	// InfiniteAreaLightIS Private Data
	Spectrum Lbase;
	MIPMap<Spectrum> *radianceMap;
	Distribution1D *uDistrib, **vDistribs;
};

// InfiniteAreaLightIS Method Definitions
InfiniteAreaLightIS::~InfiniteAreaLightIS() {
	delete radianceMap;
}
InfiniteAreaLightIS
	::InfiniteAreaLightIS(const Transform &light2world,
				const Spectrum &L, int ns,
				const string &texmap)
	: Light(light2world, ns)
{
	radianceMap = NULL;
	if (texmap != "") {
	int width, height;
	Spectrum *texels =
		ReadImage(texmap, &width, &height);
	if (texels)
		radianceMap =
		new MIPMap<Spectrum>(width, height, texels);
	// Compute scalar-valued image from environment map
	float filter = 1.f / max(width, height);
	int nu = width, nv = height;
	float *img = new float[width*height];
	for (int x = 0; x < nu; ++x) {
		float xp = (float)x / (float)nu;
		for (int y = 0; y < nv; ++y) {
			float yp = (float)y / (float)nv;
			img[y+x*nv] = radianceMap->Lookup(xp, yp, filter).y();
		}
	}
	// Initialize sampling PDFs for infinite area light
	float *func = (float *)alloca(max(nu, nv) * sizeof(float));
	float *sinVals = (float *)alloca(nv * sizeof(float));
	for (int i = 0; i < nv; ++i)
		sinVals[i] = sin(M_PI * float(i+.5)/float(nv));
	vDistribs = new Distribution1D *[nu];
	for (int u = 0; u < nu; ++u) {
		// Compute sampling distribution for column _u_
		for (int v = 0; v < nv; ++v)
			func[v] = img[u*nv+v] *= sinVals[v];
		vDistribs[u] = new Distribution1D(func, nv);
	}
	// Compute sampling distribution for columns of image
	for (int u = 0; u < nu; ++u)
		func[u] = vDistribs[u]->funcInt;
	uDistrib = new Distribution1D(func, nu);
	delete[] img;
	delete[] texels;
	}
	Lbase = L;
}
Spectrum
	InfiniteAreaLightIS::Le(const RayDifferential &r) const {
	Vector w = r.d;
	// Compute infinite light radiance for direction
	Spectrum L = Lbase;
	if (radianceMap != NULL) {
		Vector wh = Normalize(WorldToLight(w));
		float s = SphericalPhi(wh) * INV_TWOPI;
		float t = SphericalTheta(wh) * INV_PI;
		L *= radianceMap->Lookup(s, t);
	}
	return L;
}
Spectrum InfiniteAreaLightIS::Sample_L(const Point &p, float u1,
		float u2, Vector *wi, float *pdf,
		VisibilityTester *visibility) const {
	// Find floating-point $(u,v)$ sample coordinates
	float pdfs[2];
	float fu = uDistrib->Sample(u1, &pdfs[0]);
	int u = Clamp(Float2Int(fu), 0, uDistrib->count-1);
	float fv = vDistribs[u]->Sample(u2, &pdfs[1]);
	if (pdfs[0] == 0.f || pdfs[1] == 0.f) {
		*pdf = 0.f;
		return Spectrum(0.f);
	}
	// Convert sample point to direction on the unit sphere
	float theta = fv * vDistribs[u]->invCount * M_PI;
	float phi = fu * uDistrib->invCount * 2.f * M_PI;
	float costheta = cos(theta), sintheta = sin(theta);
	if (sintheta == 0.f) return 0.f;
	float sinphi = sin(phi), cosphi = cos(phi);
	*wi = LightToWorld(Vector(sintheta * cosphi, sintheta * sinphi,
	                          costheta));
	// Compute PDF for sampled direction
	*pdf = (pdfs[0] * pdfs[1]) / (2. * M_PI * M_PI * sintheta);
	// Return radiance value for direction
	visibility->SetRay(p, *wi);
	return Lbase * radianceMap->Lookup(fu * uDistrib->invCount,
		fv * vDistribs[u]->invCount);
}
float InfiniteAreaLightIS::Pdf(const Point &,
		const Vector &w) const {
	Vector wi = WorldToLight(w);
	float theta = SphericalTheta(wi), phi = SphericalPhi(wi);
	int u = Clamp(Float2Int(phi * INV_TWOPI * uDistrib->count),
                  0, uDistrib->count-1);
	int v = Clamp(Float2Int(theta * INV_PI * vDistribs[u]->count),
                  0, vDistribs[u]->count-1);
	return (uDistrib->func[u] * vDistribs[u]->func[v]) /
           (uDistrib->funcInt * vDistribs[u]->funcInt) *
           1.f / (2.f * M_PI * M_PI * sin(theta));
}
Spectrum InfiniteAreaLightIS::Sample_L(const Scene *scene,
		float u1, float u2, float u3, float u4,
		Ray *ray, float *pdf) const {
	// Choose two points _p1_ and _p2_ on scene bounding sphere
	Point worldCenter;
	float worldRadius;
	scene->WorldBound().BoundingSphere(&worldCenter,
		&worldRadius);
	worldRadius *= 1.01f;
	Point p1 = worldCenter + worldRadius *
		UniformSampleSphere(u1, u2);
	Point p2 = worldCenter + worldRadius *
		UniformSampleSphere(u3, u4);
	// Construct ray between _p1_ and _p2_
	ray->o = p1;
	ray->d = Normalize(p2-p1);
	// Compute _InfiniteAreaLightIS_ ray weight
	Vector to_center = Normalize(worldCenter - p1);
	float costheta = AbsDot(to_center,ray->d);
	*pdf =
		costheta / ((4.f * M_PI * worldRadius * worldRadius));
	return Le(RayDifferential(ray->o, -ray->d));
}
Spectrum InfiniteAreaLightIS::Sample_L(const Point &p,
		Vector *wi, VisibilityTester *visibility) const {
	float pdf;
	Spectrum L = Sample_L(p, RandomFloat(), RandomFloat(),
		wi, &pdf, visibility);
	if (pdf == 0.) return Spectrum(0.);
	return L / pdf;
}

extern "C" DLLEXPORT Light *CreateLight(const Transform &light2world,
		const ParamSet &paramSet) {
	Spectrum L = paramSet.FindOneSpectrum("L", Spectrum(1.0));
	string texmap = paramSet.FindOneString("mapname", "");
	int nSamples = paramSet.FindOneInt("nsamples", 1);

	return new InfiniteAreaLightIS(light2world, L, nSamples, texmap);
}
