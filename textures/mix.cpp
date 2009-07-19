
/*
 * pbrt source code Copyright(c) 1998-2007 Matt Pharr and Greg Humphreys
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 * (See file License.txt for complete license)
 */

// mix.cpp*
#include "pbrt.h"
#include "texture.h"
#include "paramset.h"
// MixTexture Declarations
template <class T>
class MixTexture : public Texture<T> {
public:
	// MixTexture Public Methods
	MixTexture(Reference<Texture<T> > t1,
			   Reference<Texture<T> > t2,
			   Reference<Texture<float> > amt) {
		tex1 = t1;
		tex2 = t2;
		amount = amt;
	}
	T Evaluate(const DifferentialGeometry &dg) const {
		T t1 = tex1->Evaluate(dg), t2 = tex2->Evaluate(dg);
		float amt = amount->Evaluate(dg);
		return (1.f - amt) * t1 + amt * t2;
	}
private:
	Reference<Texture<T> > tex1, tex2;
	Reference<Texture<float> > amount;
};
// MixTexture Method Definitions
extern "C" DLLEXPORT Texture<float> * CreateFloatTexture(const Transform &tex2world,
		const TextureParams &tp) {
	return new MixTexture<float>(
		tp.GetFloatTexture("tex1", 0.f),
		tp.GetFloatTexture("tex2", 1.f),
		tp.GetFloatTexture("amount", 0.5f));
}

extern "C" DLLEXPORT Texture<Spectrum> * CreateSpectrumTexture(const Transform &tex2world,
		const TextureParams &tp) {
	return new MixTexture<Spectrum>(
		tp.GetSpectrumTexture("tex1", 0.f),
		tp.GetSpectrumTexture("tex2", 1.f),
		tp.GetFloatTexture("amount", 0.5f));
}
