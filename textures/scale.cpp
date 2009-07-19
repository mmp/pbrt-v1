
/*
 * pbrt source code Copyright(c) 1998-2007 Matt Pharr and Greg Humphreys
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 * (See file License.txt for complete license)
 */

// scale.cpp*
#include "pbrt.h"
#include "texture.h"
#include "paramset.h"
// ScaleTexture Declarations
template <class T1, class T2>
class ScaleTexture : public Texture<T2> {
public:
	// ScaleTexture Public Methods
	ScaleTexture(Reference<Texture<T1> > t1,
			Reference<Texture<T2> > t2) {
		tex1 = t1;
		tex2 = t2;
	}
	T2 Evaluate(
			const DifferentialGeometry &dg) const {
		return tex1->Evaluate(dg) * tex2->Evaluate(dg);
	}
private:
	Reference<Texture<T1> > tex1;
	Reference<Texture<T2> > tex2;
};
// ScaleTexture Method Definitions
extern "C" DLLEXPORT Texture<float> * CreateFloatTexture(const Transform &tex2world,
		const TextureParams &tp) {
	return new ScaleTexture<float, float>(tp.GetFloatTexture("tex1", 1.f),
		tp.GetFloatTexture("tex2", 1.f));
}

extern "C" DLLEXPORT Texture<Spectrum> * CreateSpectrumTexture(const Transform &tex2world,
		const TextureParams &tp) {
	return new ScaleTexture<Spectrum, Spectrum>(
		tp.GetSpectrumTexture("tex1", Spectrum(1.f)),
		tp.GetSpectrumTexture("tex2", Spectrum(1.f)));
}
