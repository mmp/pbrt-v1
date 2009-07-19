
/*
 * pbrt source code Copyright(c) 1998-2007 Matt Pharr and Greg Humphreys
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 * (See file License.txt for complete license)
 */

// box.cpp*
#include "sampling.h"
#include "paramset.h"
// Box Filter Declarations
class BoxFilter : public Filter {
public:
	BoxFilter(float xw, float yw) : Filter(xw, yw) { }
	float Evaluate(float x, float y) const;
};
// Box Filter Method Definitions
float BoxFilter::Evaluate(float x, float y) const {
	return 1.;
}
extern "C" DLLEXPORT Filter *CreateFilter(const ParamSet &ps) {
	float xw = ps.FindOneFloat("xwidth", .5f);
	float yw = ps.FindOneFloat("ywidth", .5f);
	return new BoxFilter(xw, yw);
}
