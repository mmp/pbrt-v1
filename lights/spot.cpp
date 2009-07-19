
/*
 * pbrt source code Copyright(c) 1998-2007 Matt Pharr and Greg Humphreys
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 * (See file License.txt for complete license)
 */

// spot.cpp*
#include "pbrt.h"
#include "light.h"
#include "shape.h"
// SpotLight Declarations
class SpotLight : public Light {
public:
	// SpotLight Public Methods
	SpotLight(const Transform &light2world, const Spectrum &, float width, float fall);
	Spectrum Sample_L(const Point &p, Vector *wi, VisibilityTester *vis) const;
	bool IsDeltaLight() const { return true; }
	float Falloff(const Vector &w) const;
	Spectrum Power(const Scene *) const {
		return Intensity * 2.f * M_PI *
			(1.f - .5f * (cosFalloffStart + cosTotalWidth));
	}
	Spectrum Sample_L(const Point &P, float u1, float u2,
			Vector *wo, float *pdf, VisibilityTester *visibility) const;
	Spectrum Sample_L(const Scene *scene, float u1, float u2,
			float u3, float u4, Ray *ray, float *pdf) const;
	float Pdf(const Point &, const Vector &) const;
private:
	// SpotLight Private Data
	float cosTotalWidth, cosFalloffStart;
	Point lightPos;
	Spectrum Intensity;
};
// SpotLight Method Definitions
SpotLight::SpotLight(const Transform &light2world,
		const Spectrum &intensity, float width, float fall)
	: Light(light2world) {
	lightPos = LightToWorld(Point(0,0,0));
	Intensity = intensity;
	cosTotalWidth = cosf(Radians(width));
	cosFalloffStart = cosf(Radians(fall));
}
Spectrum SpotLight::Sample_L(const Point &p, Vector *wi,
		VisibilityTester *visibility) const {
	*wi = Normalize(lightPos - p);
	visibility->SetSegment(p, lightPos);
	return Intensity * Falloff(-*wi) /
		DistanceSquared(lightPos, p);
}
float SpotLight::Falloff(const Vector &w) const {
	Vector wl = Normalize(WorldToLight(w));
	float costheta = wl.z;
	if (costheta < cosTotalWidth)
		return 0.;
 	if (costheta > cosFalloffStart)
		return 1.;
	// Compute falloff inside spotlight cone
	float delta = (costheta - cosTotalWidth) /
		(cosFalloffStart - cosTotalWidth);
	return delta*delta*delta*delta;
}
Spectrum SpotLight::Sample_L(const Point &p, float u1, float u2,
		Vector *wi, float *pdf, VisibilityTester *visibility) const {
	*pdf = 1.f;
	return Sample_L(p, wi, visibility);
}
float SpotLight::Pdf(const Point &, const Vector &) const {
	return 0.;
}
Spectrum SpotLight::Sample_L(const Scene *scene, float u1,
		float u2, float u3, float u4,
		Ray *ray, float *pdf) const {
	ray->o = lightPos;
	Vector v = UniformSampleCone(u1, u2, cosTotalWidth);
	ray->d = LightToWorld(v);
	*pdf = UniformConePdf(cosTotalWidth);
	return Intensity * Falloff(ray->d);
}
extern "C" DLLEXPORT Light *CreateLight(const Transform &l2w, const ParamSet &paramSet) {
	Spectrum I = paramSet.FindOneSpectrum("I", Spectrum(1.0));
	float coneangle = paramSet.FindOneFloat("coneangle", 30.);
	float conedelta = paramSet.FindOneFloat("conedeltaangle", 5.);
	// Compute spotlight world to light transformation
	Point from = paramSet.FindOnePoint("from", Point(0,0,0));
	Point to = paramSet.FindOnePoint("to", Point(0,0,1));
	Vector dir = Normalize(to - from);
	Vector du, dv;
	CoordinateSystem(dir, &du, &dv);
	Transform dirToZ =
		Transform(new Matrix4x4( du.x,  du.y,  du.z, 0.,
	                                 dv.x,  dv.y,  dv.z, 0.,
	                                dir.x, dir.y, dir.z, 0.,
	                                    0,     0,     0, 1.));
	Transform light2world =
	l2w *
	Translate(Vector(from.x, from.y, from.z)) *
	dirToZ.GetInverse();
	return new SpotLight(light2world, I, coneangle,
		coneangle-conedelta);
}
