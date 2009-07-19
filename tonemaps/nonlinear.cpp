
/*
 * pbrt source code Copyright(c) 1998-2007 Matt Pharr and Greg Humphreys
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 * (See file License.txt for complete license)
 */

// nonlinear.cpp*
#include "tonemap.h"
#include "paramset.h"
// NonLinearOp Declarations
class NonLinearOp : public ToneMap {
public:
	// NonLinearOp Public Methods
	NonLinearOp(float my) { maxY = my; }
	void Map(const float *y, int xRes, int yRes,
			float maxDisplayY, float *scale) const {
		float invY2;
		if (maxY <= 0.f) {
			// Compute world adaptation luminance, _Ywa_
			float Ywa = 0.;
			for (int i = 0; i < xRes * yRes; ++i)
				if (y[i] > 0) Ywa += logf(y[i]);
			Ywa = expf(Ywa / (xRes * yRes));
			Ywa /= 683.f;
			invY2 = 1.f / (Ywa * Ywa);
		}
		else invY2 = 1.f / (maxY * maxY);
		for (int i = 0; i < xRes * yRes; ++i) {
			float ys = y[i] / 683.f;
			scale[i] = maxDisplayY / 683.f *
				(1.f + ys * invY2) / (1.f + ys);
		}
	}
private:
	float maxY;
};
// NonLinearOp Method Definitions
extern "C" DLLEXPORT ToneMap *CreateToneMap(const ParamSet &ps) {
	float maxy = ps.FindOneFloat("maxY", 0.f);
	return new NonLinearOp(maxy);
}
