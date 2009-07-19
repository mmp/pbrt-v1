
/*
 * pbrt source code Copyright(c) 1998-2007 Matt Pharr and Greg Humphreys
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 * (See file License.txt for complete license)
 */

// bluepaint.cpp*
#include "pbrt.h"
#include "material.h"
// BluePaint Class Declarations
class BluePaint : public Material {
public:
	BluePaint(Reference<Texture<float> > bump) : bumpMap(bump) { }
	BSDF *GetBSDF(const DifferentialGeometry &dgGeom, const DifferentialGeometry &dgShading) const;
	Reference<Texture<float> > bumpMap;
};
// BluePaint Method Definitions
BSDF *BluePaint::GetBSDF(const DifferentialGeometry &dgGeom, const DifferentialGeometry &dgShading) const {
	// Declare bluepaint coefficients
	static float diffuse[3] = {  0.3094f,    0.39667f,   0.70837f  };
	static float xy0[3] =     {  0.870567f,  0.857255f,  0.670982f };
	static float z0[3] =      {  0.803624f,  0.774290f,  0.586674f };
	static float e0[3] =      { 21.820103f, 18.597755f,  7.472717f };
	static float xy1[3] =     { -0.451218f, -0.406681f, -0.477976f };
	static float z1[3] =      {  0.023123f,  0.017625f,  0.227295f };
	static float e1[3] =      {  2.774499f,  2.581499f,  3.677653f };
	static float xy2[3] =     { -1.031545f, -1.029426f, -1.026588f };
	static float z2[3] =      {  0.706734f,  0.696530f,  0.687715f };
	static float e2[3] =      { 66.899060f, 63.767912f, 57.489181f };
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
	return new BluePaint(bumpMap);
}
