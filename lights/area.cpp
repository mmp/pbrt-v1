
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

// area.cpp*
#include "light.h"
#include "primitive.h"
// AreaLight Method Definitions
AreaLight::AreaLight(const Transform &light2world,
		const Spectrum &le, int ns,
		const Reference<Shape> &s)
	: Light(light2world, ns) {
	Lemit = le;
	if (s->CanIntersect())
		shape = s;
	else {
		// Create _ShapeSet_ for _Shape_
		Reference<Shape> shapeSet = s;
		vector<Reference<Shape> > todo, done;
		todo.push_back(shapeSet);
		while (todo.size()) {
			Reference<Shape> sh = todo.back();
			todo.pop_back();
			if (sh->CanIntersect())
				done.push_back(sh);
			else
				sh->Refine(todo);
		}
		if (done.size() == 1) shape = done[0];
		else {
			if (done.size() > 16)
				Warning("Area light geometry turned into %d shapes; "
					"may be very inefficient.", (int)done.size());
			shape = new ShapeSet(done, s->ObjectToWorld, s->reverseOrientation);
		}
	}
	area = shape->Area();
}
Spectrum AreaLight::Sample_L(const Point &p,
		const Normal &n, float u1, float u2,
		Vector *wi, float *pdf,
		VisibilityTester *visibility) const {
	Normal ns;
	Point ps = shape->Sample(p, u1, u2, &ns);
	*wi = Normalize(ps - p);
	*pdf = shape->Pdf(p, *wi);
	visibility->SetSegment(p, ps);
	return L(ps, ns, -*wi);
}
float AreaLight::Pdf(const Point &p, const Normal &N,
		const Vector &wi) const {
	return shape->Pdf(p, wi);
}
Spectrum AreaLight::Sample_L(const Point &P,
		float u1, float u2, Vector *wo, float *pdf,
		VisibilityTester *visibility) const {
	Normal Ns;
	Point Ps = shape->Sample(P, u1, u2, &Ns);
	*wo = Normalize(Ps - P);
	*pdf = shape->Pdf(P, *wo);
	visibility->SetSegment(P, Ps);
	return L(Ps, Ns, -*wo);
}
Spectrum AreaLight::Sample_L(const Scene *scene, float u1,
		float u2, float u3, float u4,
		Ray *ray, float *pdf) const {
	Normal ns;
	ray->o = shape->Sample(u1, u2, &ns);
	ray->d = UniformSampleSphere(u3, u4);
	if (Dot(ray->d, ns) < 0.) ray->d *= -1;
	*pdf = shape->Pdf(ray->o) * INV_TWOPI;
	return L(ray->o, ns, ray->d);
}
float AreaLight::Pdf(const Point &P, const Vector &w) const {
	return shape->Pdf(P, w);
}
Spectrum AreaLight::Sample_L(const Point &P, Vector *wo,
		VisibilityTester *visibility) const {
	Normal Ns;
	Point Ps = shape->Sample(P, RandomFloat(), RandomFloat(), &Ns);
	*wo = Normalize(Ps - P);
	visibility->SetSegment(P, Ps);
	float pdf = shape->Pdf(P, *wo);
	if (pdf == 0.f) return Spectrum(0.f);
	return L(P, Ns, -*wo) /	pdf;
}
extern "C" DLLEXPORT Light *CreateAreaLight(const Transform &light2world, const ParamSet &paramSet,
		const Reference<Shape> &shape) {
	Spectrum L = paramSet.FindOneSpectrum("L", Spectrum(1.0));
	int nSamples = paramSet.FindOneInt("nsamples", 1);
	return new AreaLight(light2world, L, nSamples, shape);
}
