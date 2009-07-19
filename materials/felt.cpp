
/*
 * pbrt source code Copyright(c) 1998-2007 Matt Pharr and Greg Humphreys
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 * (See file License.txt for complete license)
 */

// felt.cpp*
#include "pbrt.h"
#include "material.h"
// Felt Class Declarations
class Felt : public Material {
public:
	Felt(Reference<Texture<float> > bump) : bumpMap(bump) { }
	BSDF *GetBSDF(const DifferentialGeometry &dgGeom, const DifferentialGeometry &dgShading) const;
	Reference<Texture<float> > bumpMap;
};

// Felt Method Definitions
BSDF *Felt::GetBSDF(const DifferentialGeometry &dgGeom, const DifferentialGeometry &dgShading) const {
	// Declare felt coefficients
	static float diffuse[3] = {  0.025865f,  0.025865f,  0.025865f};
	static float xy0[3] =     { -0.304075f, -0.304075f, -0.304075f};
	static float z0[3] =      { -0.065992f, -0.065992f, -0.065992f};
	static float e0[3] =      {  3.047892f,  3.047892f,  3.047892f};
	static float xy1[3] =     { -0.749561f, -0.749561f, -0.749561f};
	static float z1[3] =      { -1.167929f, -1.167929f, -1.167929f};
	static float e1[3] =      {  6.931827f,  6.931827f,  6.931827f};
	static float xy2[3] =     {  1.004921f,  1.004921f,  1.004921f};
	static float z2[3] =      { -0.205529f, -0.205529f, -0.205529f};
	static float e2[3] =      { 94.117332f, 94.117332f, 94.117332f};
	static Spectrum xy[3] = { Spectrum(xy0), Spectrum(xy1), Spectrum(xy2) };
	static Spectrum z[3] = { Spectrum(z0), Spectrum(z1), Spectrum(z2) };
	static Spectrum e[3] = { Spectrum(e0), Spectrum(e1), Spectrum(e2) };
	// Allocate _BSDF_, possibly doing bump-mapping with _bumpMap_
	DifferentialGeometry dgs;
	if (bumpMap)
		Bump(bumpMap, dgGeom, dgShading, &dgs);
	else
		dgs = dgShading;
	BSDF *bsdf = BSDF_ALLOC(BSDF)(dgs, dgGeom.nn);
	bsdf->Add(BSDF_ALLOC(Lafortune)(Spectrum(diffuse), 3, xy, xy, z, e,
		BxDFType(BSDF_REFLECTION | BSDF_DIFFUSE)));
	return bsdf;
}
extern "C" DLLEXPORT Material * CreateMaterial(const Transform &xform,
		const TextureParams &mp) {
	Reference<Texture<float> > bumpMap = mp.GetFloatTexture("bumpmap", 0.f);
	return new Felt(bumpMap);
}
