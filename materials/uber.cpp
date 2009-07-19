
/*
 * pbrt source code Copyright(c) 1998-2007 Matt Pharr and Greg Humphreys
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 * (See file License.txt for complete license)
 */

// uber.cpp*
#include "pbrt.h"
#include "material.h"
// UberMaterial Class Declarations
class UberMaterial : public Material {
public:
	// UberMaterial Method Declarations
	UberMaterial(Reference<Texture<Spectrum> > kd,
		Reference<Texture<Spectrum> > ks,
		Reference<Texture<Spectrum> > kr,
		Reference<Texture<float> > rough,
		Reference<Texture<Spectrum> > op,
		Reference<Texture<float> > bump) {
		Kd = kd;
		Ks = ks;
		Kr = kr;
		roughness = rough;
		opacity = op;
		bumpMap = bump;
	}
	BSDF *GetBSDF(const DifferentialGeometry &dgGeom, const DifferentialGeometry &dgShading) const;
private:
	// UberMaterial Private Data
	Reference<Texture<Spectrum> > Kd, Ks, Kr, opacity;
	Reference<Texture<float> > roughness;
	Reference<Texture<float> > bumpMap;
};
// UberMaterial Method Definitions
BSDF *UberMaterial::GetBSDF(const DifferentialGeometry &dgGeom, const DifferentialGeometry &dgShading) const {
	// Allocate _BSDF_, possibly doing bump-mapping with _bumpMap_
	DifferentialGeometry dgs;
	if (bumpMap)
		Bump(bumpMap, dgGeom, dgShading, &dgs);
	else
		dgs = dgShading;
	BSDF *bsdf = BSDF_ALLOC(BSDF)(dgs, dgGeom.nn);

	Spectrum op = opacity->Evaluate(dgs).Clamp();
	if (op != Spectrum(1.)) {
	    BxDF *tr = BSDF_ALLOC(SpecularTransmission)(-op + Spectrum(1.), 1., 1.);
	    bsdf->Add(tr);
	}

	Spectrum kd = op * Kd->Evaluate(dgs).Clamp();
	if (!kd.Black()) {
	    BxDF *diff = BSDF_ALLOC(Lambertian)(kd);
	    bsdf->Add(diff);
	}

	Spectrum ks = op * Ks->Evaluate(dgs).Clamp();
	if (!ks.Black()) {
	    Fresnel *fresnel = BSDF_ALLOC(FresnelDielectric)(1.5f, 1.f);
	    float rough = roughness->Evaluate(dgs);
	    BxDF *spec = BSDF_ALLOC(Microfacet)(ks, fresnel, BSDF_ALLOC(Blinn)(1.f / rough));
	    bsdf->Add(spec);
	}

	Spectrum kr = op * Kr->Evaluate(dgs).Clamp();
	if (!kr.Black()) {
		Fresnel *fresnel = BSDF_ALLOC(FresnelDielectric)(1.5f, 1.f);
		bsdf->Add(BSDF_ALLOC(SpecularReflection)(kr, fresnel));
	}

	return bsdf;
}
// UberMaterial Dynamic Creation Routine
extern "C" DLLEXPORT Material * CreateMaterial(const Transform &xform,
		const TextureParams &mp) {
	Reference<Texture<Spectrum> > Kd = mp.GetSpectrumTexture("Kd", Spectrum(1.f));
	Reference<Texture<Spectrum> > Ks = mp.GetSpectrumTexture("Ks", Spectrum(1.f));
	Reference<Texture<Spectrum> > Kr = mp.GetSpectrumTexture("Kr", Spectrum(0.f));
	Reference<Texture<float> > roughness = mp.GetFloatTexture("roughness", .1f);
	Reference<Texture<Spectrum> > opacity = mp.GetSpectrumTexture("opacity", 1.f);
	Reference<Texture<float> > bumpMap = mp.GetFloatTexture("bumpmap", 0.f);
	return new UberMaterial(Kd, Ks, Kr, roughness, opacity, bumpMap);
}
