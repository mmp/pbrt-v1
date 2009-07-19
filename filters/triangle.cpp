
/*
 * pbrt source code Copyright(c) 1998-2007 Matt Pharr and Greg Humphreys
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 * (See file License.txt for complete license)
 */

// triangle.cpp*
#include "sampling.h"
#include "paramset.h"
// Triangle Filter Declarations
class TriangleFilter : public Filter {
public:
	TriangleFilter(float xw, float yw) : Filter(xw, yw) { }
	float Evaluate(float x, float y) const;
};
// Triangle Filter Method Definitions
float TriangleFilter::Evaluate(float x, float y) const {
	return max(0.f, xWidth - fabsf(x)) *
		max(0.f, yWidth - fabsf(y));
}
extern "C" DLLEXPORT Filter *CreateFilter(const ParamSet &ps) {
	// Find common filter parameters
	float xw = ps.FindOneFloat("xwidth", 2.);
	float yw = ps.FindOneFloat("ywidth", 2.);
	return new TriangleFilter(xw, yw);
}
