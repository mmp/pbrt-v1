
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

// glass.cpp*
#include "pbrt.h"
#include "material.h"
// Glass Class Declarations
class Glass : public Material {
public:
	// Glass Public Methods
	Glass(Reference<Texture<Spectrum> > r, Reference<Texture<Spectrum> > t,
			Reference<Texture<float> > i, Reference<Texture<float> > bump) {
		Kr = r;
		Kt = t;
		index = i;
		bumpMap = bump;
	}
	BSDF *GetBSDF(const DifferentialGeometry &dgGeom, const DifferentialGeometry &dgShading) const;
private:
	// Glass Private Data
	Reference<Texture<Spectrum> > Kr, Kt;
	Reference<Texture<float> > index;
	Reference<Texture<float> > bumpMap;
};
// Glass Method Definitions
BSDF *Glass::GetBSDF(const DifferentialGeometry &dgGeom, const DifferentialGeometry &dgShading) const {
	// Allocate _BSDF_, possibly doing bump-mapping with _bumpMap_
	DifferentialGeometry dgs;
	if (bumpMap)
		Bump(bumpMap, dgGeom, dgShading, &dgs);
	else
		dgs = dgShading;
	float ior = index->Evaluate(dgs);
	BSDF *bsdf = BSDF_ALLOC(BSDF)(dgs, dgGeom.nn, ior);
	Spectrum R = Kr->Evaluate(dgs).Clamp();
	Spectrum T = Kt->Evaluate(dgs).Clamp();
	if (!R.Black())
		bsdf->Add(BSDF_ALLOC(SpecularReflection)(R,
			BSDF_ALLOC(FresnelDielectric)(1., ior)));
	if (!T.Black())
		bsdf->Add(BSDF_ALLOC(SpecularTransmission)(T, 1., ior));
	return bsdf;
}
extern "C" DLLEXPORT Material * CreateMaterial(const Transform &xform,
		const TextureParams &mp) {
	Reference<Texture<Spectrum> > Kr = mp.GetSpectrumTexture("Kr", Spectrum(1.f));
	Reference<Texture<Spectrum> > Kt = mp.GetSpectrumTexture("Kt", Spectrum(1.f));
	Reference<Texture<float> > index = mp.GetFloatTexture("index", 1.5f);
	Reference<Texture<float> > bumpMap = mp.GetFloatTexture("bumpmap", 0.f);
	return new Glass(Kr, Kt, index, bumpMap);
}
