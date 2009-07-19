
/*
 * pbrt source code Copyright(c) 1998-2007 Matt Pharr and Greg Humphreys
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 * (See file License.txt for complete license)
 */

// maxwhite.cpp*
#include "tonemap.h"
// MaxWhiteOp Declarations
class MaxWhiteOp : public ToneMap {
public:
	// MaxWhiteOp Public Methods
	void Map(const float *y, int xRes, int yRes,
	         float maxDisplayY, float *scale) const {
		// Compute maximum luminance of all pixels
		float maxY = 0.;
		for (int i = 0; i < xRes * yRes; ++i)
			maxY = max(maxY, y[i]);
		float s = maxDisplayY / maxY;
		for (int i = 0; i < xRes * yRes; ++i)
			scale[i] = s;
	}
};
// MaxWhiteOp Method Definitions
extern "C" DLLEXPORT ToneMap *CreateToneMap(const ParamSet &ps) {
	return new MaxWhiteOp;
}
