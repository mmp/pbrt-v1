
/*
 * pbrt source code Copyright(c) 1998-2007 Matt Pharr and Greg Humphreys
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 * (See file License.txt for complete license)
 */

// contrast.cpp*
#include "tonemap.h"
#include "paramset.h"
// ContrastOp Declarations
class ContrastOp : public ToneMap {
public:
	ContrastOp(float day) { displayAdaptationY = day; }
	void Map(const float *y,
	         int xRes, int yRes,
			 float maxDisplayY, float *scale) const;
	float displayAdaptationY;
};
// ContrastOp Method Definitions
void ContrastOp::Map(const float *y, int xRes, int yRes,
		float maxDisplayY, float *scale) const {
	// Compute world adaptation luminance, _Ywa_
	float Ywa = 0.;
	for (int i = 0; i < xRes * yRes; ++i)
		if (y[i] > 0) Ywa += logf(y[i]);
	Ywa = expf(Ywa / (xRes * yRes));
	// Compute contrast-preserving scalefactor, _s_
	float s = powf((1.219f + powf(displayAdaptationY, 0.4f)) /
		(1.219f + powf(Ywa, 0.4f)), 2.5f);
	for (int i = 0; i < xRes*yRes; ++i)
		scale[i] = s;
}
extern "C" DLLEXPORT ToneMap *CreateToneMap(const ParamSet &ps) {
	float day = ps.FindOneFloat("displayadaptationY", 50.f);
	return new ContrastOp(day);
}
