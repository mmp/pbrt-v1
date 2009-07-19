
/*
 * pbrt source code Copyright(c) 1998-2007 Matt Pharr and Greg Humphreys
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 * (See file License.txt for complete license)
 */

// fbm.cpp*
#include "pbrt.h"
#include "texture.h"
#include "paramset.h"
// FBmTexture Declarations
template <class T> class FBmTexture : public Texture<T> {
public:
	// FBmTexture Public Methods
	~FBmTexture() {
		delete mapping;
	}
	FBmTexture(int oct, float roughness, TextureMapping3D *map) {
		omega = roughness;
		octaves = oct;
		mapping = map;
	}
	T Evaluate(const DifferentialGeometry &dg) const {
		Vector dpdx, dpdy;
		Point P = mapping->Map(dg, &dpdx, &dpdy);
		return FBm(P, dpdx, dpdy, omega, octaves);
	}
private:
	// FBmTexture Private Data
	int octaves;
	float omega;
	TextureMapping3D *mapping;
};
// FBmTexture Method Definitions
extern "C" DLLEXPORT Texture<float> * CreateFloatTexture(const Transform &tex2world,
		const TextureParams &tp) {
	// Initialize 3D texture mapping _map_ from _tp_
	TextureMapping3D *map = new IdentityMapping3D(tex2world);
	return new FBmTexture<float>(tp.FindInt("octaves", 8),
		tp.FindFloat("roughness", .5f), map);
}

extern "C" DLLEXPORT Texture<Spectrum> * CreateSpectrumTexture(const Transform &tex2world,
		const TextureParams &tp) {
	// Initialize 3D texture mapping _map_ from _tp_
	TextureMapping3D *map = new IdentityMapping3D(tex2world);
	return new FBmTexture<Spectrum>(tp.FindInt("octaves", 8),
		tp.FindFloat("roughness", .5f), map);
}
