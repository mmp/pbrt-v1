
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

// translucent.cpp*
#include "pbrt.h"
#include "material.h"
// Translucent Class Declarations
class Translucent : public Material {
public:
	// Translucent Public Methods
	Translucent(Reference<Texture<Spectrum> > kd, Reference<Texture<Spectrum> > ks,
			Reference<Texture<float> > rough,
			Reference<Texture<Spectrum> > refl,
			Reference<Texture<Spectrum> > trans,
			Reference<Texture<float> > bump) {
		Kd = kd;
		Ks = ks;
		roughness = rough;
		reflect = refl;
		transmit = trans;
		bumpMap = bump;
	}
	BSDF *GetBSDF(const DifferentialGeometry &dgGeom, const DifferentialGeometry &dgShading) const;
private:
	// Translucent Private Data
	Reference<Texture<Spectrum> > Kd, Ks;
	Reference<Texture<float> > roughness;
	Reference<Texture<Spectrum> > reflect, transmit;
	Reference<Texture<float> > bumpMap;
};
// Translucent Method Definitions
BSDF *Translucent::GetBSDF(const DifferentialGeometry &dgGeom, const DifferentialGeometry &dgShading) const {
	// Allocate _BSDF_, possibly doing bump-mapping with _bumpMap_
	DifferentialGeometry dgs;
	if (bumpMap)
		Bump(bumpMap, dgGeom, dgShading, &dgs);
	else
		dgs = dgShading;
	BSDF *bsdf = BSDF_ALLOC(BSDF)(dgs, dgGeom.nn);
	Spectrum r = reflect->Evaluate(dgs).Clamp();
	Spectrum t = transmit->Evaluate(dgs).Clamp();
	if (r.Black() && t.Black()) return bsdf;

	Spectrum kd = Kd->Evaluate(dgs).Clamp();
	if (!kd.Black()) {
		if (!r.Black()) bsdf->Add(BSDF_ALLOC(Lambertian)(r * kd));
		if (!t.Black()) bsdf->Add(BSDF_ALLOC(BRDFToBTDF)(BSDF_ALLOC(Lambertian)(t * kd)));
	}
	Spectrum ks = Ks->Evaluate(dgs).Clamp();
	if (!ks.Black()) {
		float rough = roughness->Evaluate(dgs);
		if (!r.Black()) {
			Fresnel *fresnel = BSDF_ALLOC(FresnelDielectric)(1.5f, 1.f);
			bsdf->Add(BSDF_ALLOC(Microfacet)(r * ks, fresnel,
				BSDF_ALLOC(Blinn)(1.f / rough)));
		}
		if (!t.Black()) {
			Fresnel *fresnel = BSDF_ALLOC(FresnelDielectric)(1.5f, 1.f);
			bsdf->Add(BSDF_ALLOC(BRDFToBTDF)(BSDF_ALLOC(Microfacet)(t * ks, fresnel,
				BSDF_ALLOC(Blinn)(1.f / rough))));
		}
	}
	return bsdf;
}
extern "C" DLLEXPORT Material * CreateMaterial(const Transform &xform,
		const TextureParams &mp) {
	Reference<Texture<Spectrum> > Kd = mp.GetSpectrumTexture("Kd", Spectrum(1.f));
	Reference<Texture<Spectrum> > Ks = mp.GetSpectrumTexture("Ks", Spectrum(1.f));
	Reference<Texture<Spectrum> > reflect = mp.GetSpectrumTexture("reflect", Spectrum(0.5f));
	Reference<Texture<Spectrum> > transmit = mp.GetSpectrumTexture("transmit", Spectrum(0.5f));
	Reference<Texture<float> > roughness = mp.GetFloatTexture("roughness", .1f);
	Reference<Texture<float> > bumpMap = mp.GetFloatTexture("bumpmap", 0.f);
	return new Translucent(Kd, Ks, roughness, reflect, transmit, bumpMap);
}
