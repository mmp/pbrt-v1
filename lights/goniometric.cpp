
/*
 * pbrt source code Copyright(c) 1998-2007 Matt Pharr and Greg Humphreys
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 * (See file License.txt for complete license)
 */

// goniometric.cpp*
#include "pbrt.h"
#include "light.h"
#include "shape.h"
#include "scene.h"
#include "mipmap.h"
// GonioPhotometricLight Declarations
class GonioPhotometricLight : public Light {
public:
	// GonioPhotometricLight Public Methods
	GonioPhotometricLight(const Transform &light2world, const Spectrum &, const
	string &texname);
	Spectrum Sample_L(const Point &p, Vector *wi, VisibilityTester *vis) const;
	~GonioPhotometricLight() { delete mipmap; }
	bool IsDeltaLight() const { return true; }
	Spectrum Scale(const Vector &w) const {
		Vector wp = Normalize(WorldToLight(w));
		swap(wp.y, wp.z);
		float theta = SphericalTheta(wp);
		float phi   = SphericalPhi(wp);
		float s = phi * INV_TWOPI, t = theta * INV_PI;
		return mipmap ? mipmap->Lookup(s, t) : 1.f;
	}
	Spectrum Power(const Scene *) const {
		return 4.f * M_PI * Intensity *
			mipmap->Lookup(.5f, .5f, .5f);
	}
	Spectrum Sample_L(const Point &P, float u1, float u2, Vector *wo,
		float *pdf, VisibilityTester *visibility) const;
	Spectrum Sample_L(const Scene *scene, float u1, float u2,
			float u3, float u4, Ray *ray, float *pdf) const;
	float Pdf(const Point &, const Vector &) const;
private:
	// GonioPhotometricLight Private Data
	Point lightPos;
	Spectrum Intensity;
	MIPMap<Spectrum> *mipmap;

};
// GonioPhotometricLight Method Definitions
Spectrum GonioPhotometricLight::Sample_L(const Point &p, Vector *wi,
		VisibilityTester *visibility) const {
	*wi = Normalize(lightPos - p);
	visibility->SetSegment(p, lightPos);
	return Intensity * Scale(-*wi) / DistanceSquared(lightPos, p);
}
GonioPhotometricLight::GonioPhotometricLight(
		const Transform &light2world,
		const Spectrum &intensity, const string &texname)
	: Light(light2world) {
	lightPos = LightToWorld(Point(0,0,0));
	Intensity = intensity;
	// Create _mipmap_ for _GonioPhotometricLight_
	int width, height;
	Spectrum *texels = ReadImage(texname, &width, &height);
	if (texels) {
		mipmap = new MIPMap<Spectrum>(width, height, texels);
		delete[] texels;
	}
	else mipmap = NULL;
}
Spectrum GonioPhotometricLight::Sample_L(const Point &P, float u1, float u2,
		Vector *wo, float *pdf,
		VisibilityTester *visibility) const {
	*wo = Normalize(lightPos - P);
	*pdf = 1.f;
	visibility->SetSegment(P, lightPos);
	return Intensity * Scale(-*wo) / DistanceSquared(lightPos, P);
}
Spectrum GonioPhotometricLight::Sample_L(const Scene *scene, float u1, float u2,
		float u3, float u4, Ray *ray, float *pdf) const {
	ray->o = lightPos;
	ray->d = UniformSampleSphere(u1, u2);
	*pdf = UniformSpherePdf();
	return Intensity * Scale(ray->d);
}
float GonioPhotometricLight::Pdf(const Point &, const Vector &) const {
	return 0.;
}
extern "C" DLLEXPORT Light *CreateLight(const Transform &light2world,
		const ParamSet &paramSet) {
	Spectrum I = paramSet.FindOneSpectrum("I", Spectrum(1.0));
	string texname = paramSet.FindOneString("mapname", "");
	return new GonioPhotometricLight(light2world, I, texname);
}
