
/*
 * pbrt source code Copyright(c) 1998-2007 Matt Pharr and Greg Humphreys
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 * (See file License.txt for complete license)
 */

// mirror.cpp*
#include "pbrt.h"
#include "material.h"
// Mirror Class Declarations
class Mirror : public Material {
public:
	// Mirror Public Methods
	Mirror(Reference<Texture<Spectrum> > r, Reference<Texture<float> > bump) {
		Kr = r;
		bumpMap = bump;
	}
	BSDF *GetBSDF(const DifferentialGeometry &dgGeom, const DifferentialGeometry &dgShading) const;
private:
	// Mirror Private Data
	Reference<Texture<Spectrum> > Kr;
	Reference<Texture<float> > bumpMap;
};
// Mirror Method Definitions
BSDF *Mirror::GetBSDF(const DifferentialGeometry &dgGeom, const DifferentialGeometry &dgShading) const {
	// Allocate _BSDF_, possibly doing bump-mapping with _bumpMap_
	DifferentialGeometry dgs;
	if (bumpMap)
		Bump(bumpMap, dgGeom, dgShading, &dgs);
	else
		dgs = dgShading;
	BSDF *bsdf = BSDF_ALLOC(BSDF)(dgs, dgGeom.nn);
	Spectrum R = Kr->Evaluate(dgs).Clamp();
	if (!R.Black())
		bsdf->Add(BSDF_ALLOC(SpecularReflection)(R,
			BSDF_ALLOC(FresnelNoOp)()));
	return bsdf;
}
extern "C" DLLEXPORT Material * CreateMaterial(const Transform &xform,
		const TextureParams &mp) {
	Reference<Texture<Spectrum> > Kr = mp.GetSpectrumTexture("Kr", Spectrum(1.f));
	Reference<Texture<float> > bumpMap = mp.GetFloatTexture("bumpmap", 0.f);
	return new Mirror(Kr, bumpMap);
}
