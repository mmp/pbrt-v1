
/*
 * pbrt source code Copyright(c) 1998-2007 Matt Pharr and Greg Humphreys
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 * (See file License.txt for complete license)
 */

// constant.cpp*
#include "pbrt.h"
#include "texture.h"
#include "paramset.h"
// ConstantTexture Method Definitions
extern "C" DLLEXPORT Texture<float> * CreateFloatTexture(const Transform &tex2world,
		const TextureParams &tp) {
	return new ConstantTexture<float>(tp.FindFloat("value", 1.f));
}

extern "C" DLLEXPORT Texture<Spectrum> * CreateSpectrumTexture(const Transform &tex2world,
		const TextureParams &tp) {
	return new ConstantTexture<Spectrum>(tp.FindSpectrum("value", Spectrum(1.f)));
}
