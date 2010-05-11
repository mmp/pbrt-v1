
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

// clay.cpp*
#include "pbrt.h"
#include "material.h"
// Clay Class Declarations
class Clay : public Material {
public:
	Clay(Reference<Texture<float> > bump) : bumpMap(bump) { }
	BSDF *GetBSDF(const DifferentialGeometry &dgGeom, const DifferentialGeometry &dgShading) const; // NOBOOK
private:
	Reference<Texture<float> > bumpMap;
};

// Clay Method Definitions
BSDF *Clay::GetBSDF(const DifferentialGeometry &dgGeom, const DifferentialGeometry &dgShading) const {
	// Declare clay coefficients
	static float diffuse[3] = {   0.383626f,   0.260749f,   0.274207f };
	static float xy0[3] =     {  -1.089701f,  -1.102701f,  -1.107603f };
	static float z0[3] =      {  -1.354682f,  -2.714801f,  -1.569866f };
	static float e0[3] =      {  17.968505f,  11.024489f,  21.270282f };
	static float xy1[3] =     {  -0.733381f,  -0.793320f,  -0.848206f };
	static float z1[3] =      {   0.676108f,   0.679314f,   0.726031f };
	static float e1[3] =      {   8.219745f,   9.055139f,  11.261951f };
	static float xy2[3] =     {  -1.010548f,  -1.012378f,  -1.011263f };
	static float z2[3] =      {   0.910783f,   0.885239f,   0.892451f };
	static float e2[3] =      { 152.912795f, 141.937171f, 201.046802f };
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
	return new Clay(bumpMap);
}
