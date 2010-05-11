
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

// primer.cpp*
#include "pbrt.h"
#include "material.h"
// Primer Class Declarations
class Primer : public Material {
public:
	Primer(Reference<Texture<float> > bump) : bumpMap(bump) { }
	BSDF *GetBSDF(const DifferentialGeometry &dgGeom, const DifferentialGeometry &dgShading) const;
	Reference<Texture<float> > bumpMap;
};

// Primer Method Definitions
BSDF *Primer::GetBSDF(const DifferentialGeometry &dgGeom, const DifferentialGeometry &dgShading) const {
	// Declare primer coefficients
	static float diffuse[3] = {  0.118230f,  0.121218f,  0.133209f};
	static float xy0[3] =     { -0.399286f, -1.033473f, -1.058104f};
	static float z0[3] =      {  0.167504f,  0.009545f, -0.068002f};
	static float e0[3] =      {  2.466633f,  7.637253f,  8.117645f};
	static float xy1[3] =     { -1.041861f, -1.100108f, -1.087779f};
	static float z1[3] =      {  0.014375f, -0.198147f, -0.053605f};
	static float e1[3] =      {  7.993722f, 29.446268f, 41.988990f};
	static float xy2[3] =     { -1.098605f, -0.379883f, -0.449038f};
	static float z2[3] =      { -0.145110f,  0.159127f,  0.173224f};
	static float e2[3] =      { 31.899719f,  2.372852f,  2.636161f};
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
	return new Primer(bumpMap);
}
