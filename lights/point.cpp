
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

// point.cpp*
#include "pbrt.h"
#include "light.h"
#include "shape.h"
// PointLight Classes
class PointLight : public Light {
public:
	// PointLight Public Methods
	PointLight(const Transform &light2world, const Spectrum &intensity);
	Spectrum Sample_L(const Point &p, Vector *wi, VisibilityTester *vis) const;
	Spectrum Power(const Scene *) const {
		return Intensity * 4.f * M_PI;
	}
	bool IsDeltaLight() const { return true; }
	Spectrum Sample_L(const Point &P, float u1, float u2,
			Vector *wo, float *pdf, VisibilityTester *visibility) const;
	Spectrum Sample_L(const Scene *scene, float u1, float u2,
			float u3, float u4, Ray *ray, float *pdf) const;
	float Pdf(const Point &, const Vector &) const;
private:
	// PointLight Private Data
	Point lightPos;
	Spectrum Intensity;
};
// PointLight Method Definitions
PointLight::PointLight(const Transform &light2world,
		const Spectrum &intensity)
	: Light(light2world) {
	lightPos = LightToWorld(Point(0,0,0));
	Intensity = intensity;
}
Spectrum PointLight::Sample_L(const Point &p, Vector *wi,
		VisibilityTester *visibility) const {
	*wi = Normalize(lightPos - p);
	visibility->SetSegment(p, lightPos);
	return Intensity / DistanceSquared(lightPos, p);
}
Spectrum PointLight::Sample_L(const Point &p, float u1,
		float u2, Vector *wi, float *pdf,
		VisibilityTester *visibility) const {
	*pdf = 1.f;
	return Sample_L(p, wi, visibility);
}
float PointLight::Pdf(const Point &, const Vector &) const {
	return 0.;
}
Spectrum PointLight::Sample_L(const Scene *scene, float u1,
		float u2, float u3, float u4,
		Ray *ray, float *pdf) const {
	ray->o = lightPos;
	ray->d = UniformSampleSphere(u1, u2);
	*pdf = UniformSpherePdf();
	return Intensity;
}
extern "C" DLLEXPORT Light *CreateLight(const Transform &light2world,
		const ParamSet &paramSet) {
	Spectrum I = paramSet.FindOneSpectrum("I", Spectrum(1.0));
	Point P = paramSet.FindOnePoint("from", Point(0,0,0));
	Transform l2w = Translate(Vector(P.x, P.y, P.z)) * light2world;
	return new PointLight(l2w, I);
}
