
/*
 * pbrt source code Copyright(c) 1998-2007 Matt Pharr and Greg Humphreys
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 * (See file License.txt for complete license)
 */

// plastic.cpp*
#include "pbrt.h"
#include "material.h"
// Plastic Class Declarations
class Plastic : public Material {
public:
	// Plastic Public Methods
	Plastic(Reference<Texture<Spectrum> > kd,
			Reference<Texture<Spectrum> > ks,
			Reference<Texture<float> > rough,
			Reference<Texture<float> > bump) {
		Kd = kd;
		Ks = ks;
		roughness = rough;
		bumpMap = bump;
	}
	BSDF *GetBSDF(const DifferentialGeometry &dgGeom,
	              const DifferentialGeometry &dgShading) const;
private:
	// Plastic Private Data
	Reference<Texture<Spectrum> > Kd, Ks;
	Reference<Texture<float> > roughness, bumpMap;
};
// Plastic Method Definitions
BSDF *Plastic::GetBSDF(const DifferentialGeometry &dgGeom,
		const DifferentialGeometry &dgShading) const {
	// Allocate _BSDF_, possibly doing bump-mapping with _bumpMap_
	DifferentialGeometry dgs;
	if (bumpMap)
		Bump(bumpMap, dgGeom, dgShading, &dgs);
	else
		dgs = dgShading;
	BSDF *bsdf = BSDF_ALLOC(BSDF)(dgs, dgGeom.nn);
	Spectrum kd = Kd->Evaluate(dgs).Clamp();
	BxDF *diff = BSDF_ALLOC(Lambertian)(kd);
	Fresnel *fresnel =
		BSDF_ALLOC(FresnelDielectric)(1.5f, 1.f);
	Spectrum ks = Ks->Evaluate(dgs).Clamp();
	float rough = roughness->Evaluate(dgs);
	BxDF *spec = BSDF_ALLOC(Microfacet)(ks, fresnel,
		BSDF_ALLOC(Blinn)(1.f / rough));
	bsdf->Add(diff);
	bsdf->Add(spec);
	return bsdf;
}
// Plastic Dynamic Creation Routine
extern "C" DLLEXPORT Material * CreateMaterial(const Transform &xform,
		const TextureParams &mp) {
	Reference<Texture<Spectrum> > Kd = mp.GetSpectrumTexture("Kd", Spectrum(1.f));
	Reference<Texture<Spectrum> > Ks = mp.GetSpectrumTexture("Ks", Spectrum(1.f));
	Reference<Texture<float> > roughness = mp.GetFloatTexture("roughness", .1f);
	Reference<Texture<float> > bumpMap = mp.GetFloatTexture("bumpmap", 0.f);
	return new Plastic(Kd, Ks, roughness, bumpMap);
}
