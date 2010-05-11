
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

// substrate.cpp*
#include "pbrt.h"
#include "material.h"
// Substrate Class Declarations
class Substrate : public Material {
public:
	// Substrate Public Methods
	Substrate(Reference<Texture<Spectrum> > kd, Reference<Texture<Spectrum> > ks,
			Reference<Texture<float> > u, Reference<Texture<float> > v,
			Reference<Texture<float> > bump) {
		Kd = kd;
		Ks = ks;
		nu = u;
		nv = v;
		bumpMap = bump;
	}
	BSDF *GetBSDF(const DifferentialGeometry &dgGeom, const DifferentialGeometry &dgShading) const;
private:
	// Substrate Private Data
	Reference<Texture<Spectrum> > Kd, Ks;
	Reference<Texture<float> > nu, nv;
	Reference<Texture<float> > bumpMap;
};
// Substrate Method Definitions
BSDF *Substrate::GetBSDF(const DifferentialGeometry &dgGeom, const DifferentialGeometry &dgShading) const {
	// Allocate _BSDF_, possibly doing bump-mapping with _bumpMap_
	DifferentialGeometry dgs;
	if (bumpMap)
		Bump(bumpMap, dgGeom, dgShading, &dgs);
	else
		dgs = dgShading;
	BSDF *bsdf = BSDF_ALLOC(BSDF)(dgs, dgGeom.nn);
	Spectrum d = Kd->Evaluate(dgs).Clamp();
	Spectrum s = Ks->Evaluate(dgs).Clamp();
	float u = nu->Evaluate(dgs);
	float v = nv->Evaluate(dgs);

	bsdf->Add(BSDF_ALLOC(FresnelBlend)(d, s, BSDF_ALLOC(Anisotropic)(1.f/u, 1.f/v)));
	return bsdf;
}
extern "C" DLLEXPORT Material * CreateMaterial(const Transform &xform,
		const TextureParams &mp) {
	Reference<Texture<Spectrum> > Kd = mp.GetSpectrumTexture("Kd", Spectrum(.5f));
	Reference<Texture<Spectrum> > Ks = mp.GetSpectrumTexture("Ks", Spectrum(.5f));
	Reference<Texture<float> > uroughness = mp.GetFloatTexture("uroughness", .1f);
	Reference<Texture<float> > vroughness = mp.GetFloatTexture("vroughness", .1f);
	Reference<Texture<float> > bumpMap = mp.GetFloatTexture("bumpmap", 0.f);
	return new Substrate(Kd, Ks, uroughness, vroughness, bumpMap);
}
