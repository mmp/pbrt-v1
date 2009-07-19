
/*
 * pbrt source code Copyright(c) 1998-2007 Matt Pharr and Greg Humphreys
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 * (See file License.txt for complete license)
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
	BSDF *bsdf = BSDF_ALLOC(BSDF)(dgs, dgGeom.nn);
	Spectrum R = Kr->Evaluate(dgs).Clamp();
	Spectrum T = Kt->Evaluate(dgs).Clamp();
	float ior = index->Evaluate(dgs);
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
