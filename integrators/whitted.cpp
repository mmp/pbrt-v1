
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

// whitted.cpp*
#include "pbrt.h"
#include "transport.h"
#include "scene.h"
// WhittedIntegrator Declarations
class WhittedIntegrator : public SurfaceIntegrator {
public:
	// WhittedIntegrator Public Methods
	Spectrum Li(const Scene *scene, const RayDifferential &ray,
			const Sample *sample, float *alpha) const;
	WhittedIntegrator(int md) {
		maxDepth = md;
		rayDepth = 0;
	}
private:
	// WhittedIntegrator Private Data
	int maxDepth;
	mutable int rayDepth;
};
// WhittedIntegrator Method Definitions
Spectrum WhittedIntegrator::Li(const Scene *scene,
		const RayDifferential &ray, const Sample *sample,
		float *alpha) const {
	Intersection isect;
	Spectrum L(0.);
	bool hitSomething;
	// Search for ray-primitive intersection
	hitSomething = scene->Intersect(ray, &isect);
	if (!hitSomething) {
		// Handle ray with no intersection
		if (alpha) *alpha = 0.;
		for (u_int i = 0; i < scene->lights.size(); ++i)
			L += scene->lights[i]->Le(ray);
		if (alpha && !L.Black()) *alpha = 1.;
		return L;
	}
	else {
		// Initialize _alpha_ for ray hit
		if (alpha) *alpha = 1.;
		// Compute emitted and reflected light at ray intersection point
		// Evaluate BSDF at hit point
		BSDF *bsdf = isect.GetBSDF(ray);
		// Initialize common variables for Whitted integrator
		const Point &p = bsdf->dgShading.p;
		const Normal &n = bsdf->dgShading.nn;
		Vector wo = -ray.d;
		// Compute emitted light if ray hit an area light source
		L += isect.Le(wo);
		// Add contribution of each light source
		Vector wi;
		for (u_int i = 0; i < scene->lights.size(); ++i) {
			VisibilityTester visibility;
			Spectrum Li = scene->lights[i]->Sample_L(p, &wi, &visibility);
			if (Li.Black()) continue;
			Spectrum f = bsdf->f(wo, wi);
			if (!f.Black() && visibility.Unoccluded(scene))
				L += f * Li * AbsDot(wi, n) * visibility.Transmittance(scene);
		}
		if (rayDepth++ < maxDepth) {
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
		--rayDepth;
	}
	return L;
}
extern "C" DLLEXPORT SurfaceIntegrator *CreateSurfaceIntegrator(const ParamSet &params)
{
	int maxDepth = params.FindOneInt("maxdepth", 5);
	return new WhittedIntegrator(maxDepth);
}
