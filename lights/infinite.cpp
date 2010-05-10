
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

// infinite.cpp*
#include "pbrt.h"
#include "light.h"
#include "texture.h"
#include "shape.h"
#include "scene.h"
#include "mipmap.h"
// InfiniteAreaLight Declarations
class InfiniteAreaLight : public Light {
public:
	// InfiniteAreaLight Public Methods
	InfiniteAreaLight(const Transform &light2world,	const Spectrum &power, int ns, const string &texmap);
	~InfiniteAreaLight();
	Spectrum Power(const Scene *scene) const {
		Point worldCenter;
		float worldRadius;
		scene->WorldBound().BoundingSphere(&worldCenter,
		                                    &worldRadius);
		return Lbase * (radianceMap ? radianceMap->Lookup(.5f, .5f, .5f) : Spectrum(1.f)) *
			M_PI * worldRadius * worldRadius;
	}
	bool IsDeltaLight() const { return false; }
	Spectrum Le(const RayDifferential &r) const;
	Spectrum Sample_L(const Point &p, const Normal &n,
		float u1, float u2, Vector *wi, float *pdf,
		VisibilityTester *visibility) const;
	Spectrum Sample_L(const Point &p, float u1, float u2, Vector *wi, float *pdf,
		VisibilityTester *visibility) const;
	Spectrum Sample_L(const Scene *scene, float u1, float u2,
			float u3, float u4, Ray *ray, float *pdf) const;
	float Pdf(const Point &, const Normal &, const Vector &) const;
	float Pdf(const Point &, const Vector &) const;
	Spectrum Sample_L(const Point &P, Vector *w, VisibilityTester *visibility) const;
private:
	// InfiniteAreaLight Private Data
	Spectrum Lbase;
	MIPMap<Spectrum> *radianceMap;
};
// InfiniteAreaLight Method Definitions
InfiniteAreaLight::~InfiniteAreaLight() {
	delete radianceMap;
}
InfiniteAreaLight
	::InfiniteAreaLight(const Transform &light2world,
		                const Spectrum &L, int ns,
						const string &texmap)
	: Light(light2world, ns) {
	radianceMap = NULL;
	if (texmap != "") {
		int width, height;
		Spectrum *texels =
			ReadImage(texmap, &width, &height);
		if (texels)
			radianceMap =
				new MIPMap<Spectrum>(width, height, texels);
		delete[] texels;
	}
	Lbase = L;
}
Spectrum
	InfiniteAreaLight::Le(const RayDifferential &r) const {
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
Spectrum InfiniteAreaLight::Sample_L(const Point &p,
		const Normal &n, float u1, float u2,
		Vector *wi, float *pdf,
		VisibilityTester *visibility) const {
	// Sample cosine-weighted direction on unit sphere
	float x, y, z;
	ConcentricSampleDisk(u1, u2, &x, &y);
	z = sqrtf(max(0.f, 1.f - x*x - y*y));
	if (RandomFloat() < .5) z *= -1;
	*wi = Vector(x, y, z);
	// Compute _pdf_ for cosine-weighted infinite light direction
	*pdf = fabsf(wi->z) * INV_TWOPI;
	// Transform direction to world space
	Vector v1, v2;
	CoordinateSystem(Normalize(Vector(n)), &v1, &v2);
	*wi = Vector(v1.x * wi->x + v2.x * wi->y + n.x * wi->z,
	             v1.y * wi->x + v2.y * wi->y + n.y * wi->z,
	             v1.z * wi->x + v2.z * wi->y + n.z * wi->z);
	visibility->SetRay(p, *wi);
	return Le(RayDifferential(p, *wi));
}
float InfiniteAreaLight::Pdf(const Point &, const Normal &n,
		const Vector &wi) const {
	return AbsDot(n, wi) * INV_TWOPI;
}
Spectrum InfiniteAreaLight::Sample_L(const Point &p,
		float u1, float u2, Vector *wi, float *pdf,
		VisibilityTester *visibility) const {
	*wi = UniformSampleSphere(u1, u2);
	*pdf = UniformSpherePdf();
	visibility->SetRay(p, *wi);
	return Le(RayDifferential(p, *wi));
}
float InfiniteAreaLight::Pdf(const Point &, const Vector &) const {
	return 1.f / (4.f * M_PI);
}
Spectrum InfiniteAreaLight::Sample_L(const Scene *scene,
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
	// Compute _InfiniteAreaLight_ ray weight
	Vector to_center = Normalize(worldCenter - p1);
	float costheta = AbsDot(to_center,ray->d);
	*pdf =
		costheta / ((4.f * M_PI * worldRadius * worldRadius));
	return Le(RayDifferential(ray->o, -ray->d));
}
Spectrum InfiniteAreaLight::Sample_L(const Point &p,
		Vector *wi, VisibilityTester *visibility) const {
	float pdf;
	Spectrum L = Sample_L(p, RandomFloat(), RandomFloat(),
		wi, &pdf, visibility);
	if (pdf == 0.f) return Spectrum(0.f);
	return L / pdf;
}
extern "C" DLLEXPORT Light *CreateLight(const Transform &light2world,
		const ParamSet &paramSet) {
	Spectrum L = paramSet.FindOneSpectrum("L", Spectrum(1.0));
	string texmap = paramSet.FindOneString("mapname", "");
	int nSamples = paramSet.FindOneInt("nsamples", 1);
	return new InfiniteAreaLight(light2world, L, nSamples, texmap);
}
