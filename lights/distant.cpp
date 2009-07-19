
/*
 * pbrt source code Copyright(c) 1998-2007 Matt Pharr and Greg Humphreys
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 * (See file License.txt for complete license)
 */

// distant.cpp*
#include "pbrt.h"
#include "light.h"
#include "shape.h"
#include "scene.h"
// DistantLight Declarations
class DistantLight : public Light {
public:
	// DistantLight Public Methods
	DistantLight(const Transform &light2world, const Spectrum &radiance, const Vector &dir);
	bool IsDeltaLight() const { return true; }
	Spectrum Sample_L(const Point &p, Vector *wi, VisibilityTester *) const;
	Spectrum Power(const Scene *scene) const {
		Point worldCenter;
		float worldRadius;
		scene->WorldBound().BoundingSphere(&worldCenter,
		                                   &worldRadius);
		return L * M_PI * worldRadius * worldRadius;
	}
	Spectrum Sample_L(const Point &P, float u1, float u2, Vector *wo, float *pdf,
		VisibilityTester *visibility) const;
	Spectrum Sample_L(const Scene *scene, float u1, float u2,
		float u3, float u4, Ray *ray, float *pdf) const;
	float Pdf(const Point &, const Vector &) const;
private:
	// DistantLight Private Data
	Vector lightDir;
	Spectrum L;
};
// DistantLight Method Definitions
DistantLight::DistantLight(const Transform &light2world,
		const Spectrum &radiance, const Vector &dir)
	: Light(light2world) {
	lightDir = Normalize(LightToWorld(dir));
	L = radiance;
}
Spectrum DistantLight::Sample_L(const Point &p,
		Vector *wi, VisibilityTester *visibility) const {
	*wi = lightDir;
	visibility->SetRay(p, *wi);
	return L;
}
Spectrum DistantLight::Sample_L(const Point &p, float u1, float u2,
		Vector *wi, float *pdf, VisibilityTester *visibility) const {
	*pdf = 1.f;
	return Sample_L(p, wi, visibility);
}
float DistantLight::Pdf(const Point &, const Vector &) const {
	return 0.;
}
Spectrum DistantLight::Sample_L(const Scene *scene,
		float u1, float u2, float u3, float u4,
		Ray *ray, float *pdf) const {
	// Choose point on disk oriented toward infinite light direction
	Point worldCenter;
	float worldRadius;
	scene->WorldBound().BoundingSphere(&worldCenter,
	                                   &worldRadius);
	Vector v1, v2;
	CoordinateSystem(lightDir, &v1, &v2);
	float d1, d2;
	ConcentricSampleDisk(u1, u2, &d1, &d2);
	Point Pdisk =
		worldCenter + worldRadius * (d1 * v1 + d2 * v2);
	// Set ray origin and direction for infinite light ray
	ray->o = Pdisk + worldRadius * lightDir;
	ray->d = -lightDir;
	*pdf = 1.f / (M_PI * worldRadius * worldRadius);
	return L;
}
extern "C" DLLEXPORT Light *CreateLight(const Transform &light2world,
		const ParamSet &paramSet) {
	Spectrum L = paramSet.FindOneSpectrum("L", Spectrum(1.0));
	Point from = paramSet.FindOnePoint("from", Point(0,0,0));
	Point to = paramSet.FindOnePoint("to", Point(0,0,1));
	Vector dir = from-to;
	return new DistantLight(light2world, L, dir);
}
