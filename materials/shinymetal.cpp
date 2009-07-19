
/*
 * pbrt source code Copyright(c) 1998-2007 Matt Pharr and Greg Humphreys
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 * (See file License.txt for complete license)
 */

// shinymetal.cpp*
#include "pbrt.h"
#include "material.h"
// ShinyMetal Class Declarations
class ShinyMetal : public Material {
public:
	// ShinyMetal Public Methods
	ShinyMetal(Reference<Texture<Spectrum> > ks, Reference<Texture<float> > rough,
			Reference<Texture<Spectrum> > kr, Reference<Texture<float> > bump) {
		Ks = ks;
		roughness = rough;
		Kr = kr;
		bumpMap = bump;
	}
	BSDF *GetBSDF(const DifferentialGeometry &dgGeom, const DifferentialGeometry &dgShading) const;
private:
	// ShinyMetal Private Data
	Reference<Texture<Spectrum> > Ks, Kr;
	Reference<Texture<float> > roughness;
	Reference<Texture<float> > bumpMap;
};
// ShinyMetal Method Definitions
BSDF *ShinyMetal::GetBSDF(const DifferentialGeometry &dgGeom, const DifferentialGeometry &dgShading) const {
	// Allocate _BSDF_, possibly doing bump-mapping with _bumpMap_
	DifferentialGeometry dgs;
	if (bumpMap)
		Bump(bumpMap, dgGeom, dgShading, &dgs);
	else
		dgs = dgShading;
	BSDF *bsdf = BSDF_ALLOC(BSDF)(dgs, dgGeom.nn);
	Spectrum spec = Ks->Evaluate(dgs).Clamp();
	float rough = roughness->Evaluate(dgs);
	Spectrum R = Kr->Evaluate(dgs).Clamp();

	MicrofacetDistribution *md = BSDF_ALLOC(Blinn)(1.f / rough);
	Spectrum k = 0.;
	Fresnel *frMf = BSDF_ALLOC(FresnelConductor)(FresnelApproxEta(spec), k);
	Fresnel *frSr = BSDF_ALLOC(FresnelConductor)(FresnelApproxEta(R), k);
	bsdf->Add(BSDF_ALLOC(Microfacet)(1., frMf, md));
	bsdf->Add(BSDF_ALLOC(SpecularReflection)(1., frSr));
	return bsdf;
}
extern "C" DLLEXPORT Material * CreateMaterial(const Transform &xform,
		const TextureParams &mp) {
	Reference<Texture<Spectrum> > Kr = mp.GetSpectrumTexture("Kr", Spectrum(1.f));
	Reference<Texture<Spectrum> > Ks = mp.GetSpectrumTexture("Ks", Spectrum(1.f));
	Reference<Texture<float> > roughness = mp.GetFloatTexture("roughness", .1f);
	Reference<Texture<float> > bumpMap = mp.GetFloatTexture("bumpmap", 0.f);
	return new ShinyMetal(Ks, roughness, Kr, bumpMap);
}
