
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

// exrio.cpp*
#ifdef WIN32
#define hypotf hypot // For the OpenEXR headers
#endif

#include <ImfInputFile.h>
#include <ImfOutputFile.h>
#include <ImfChannelList.h>
#include <ImfFrameBuffer.h>
#include <ImfRgbaFile.h>
#include <half.h>
#include "pbrt.h"
#include "color.h"
using namespace Imf;
using namespace Imath;
// EXR Function Definitions
COREDLL Spectrum *ReadImage(const string &name, int *width, int *height) {
	try {
	InputFile file(name.c_str());
	Box2i dw = file.header().dataWindow();
	*width  = dw.max.x - dw.min.x + 1;
	*height = dw.max.y - dw.min.y + 1;

	half *rgb = new half[3 * *width * *height];

	FrameBuffer frameBuffer;
	frameBuffer.insert("R", Slice(HALF, (char *)rgb,
		3*sizeof(half), *width * 3 * sizeof(half), 1, 1, 0.0));
	frameBuffer.insert("G", Slice(HALF, (char *)rgb+sizeof(half),
		3*sizeof(half), *width * 3 * sizeof(half), 1, 1, 0.0));
	frameBuffer.insert("B", Slice(HALF, (char *)rgb+2*sizeof(half),
		3*sizeof(half), *width * 3 * sizeof(half), 1, 1, 0.0));

	file.setFrameBuffer(frameBuffer);
	file.readPixels(dw.min.y, dw.max.y);

	Spectrum *ret = new Spectrum[*width * *height];
	// XXX should do real RGB -> Spectrum conversion here
	for (int i = 0; i < *width * *height; ++i) {
		float c[3] = { rgb[3*i], rgb[3*i+1], rgb[3*i+2] };
		ret[i] = Spectrum(c);
	}
	delete[] rgb;
	return ret;
	} catch (const std::exception &e) {
		Error("Unable to read image file \"%s\": %s", name.c_str(),
			e.what());
		return NULL;
	}
}

COREDLL void WriteRGBAImage(const string &name, float *pixels,
		float *alpha, int xRes, int yRes,
		int totalXRes, int totalYRes,
		int xOffset, int yOffset) {
    Rgba *hrgba = new Rgba[xRes * yRes];
    for (int i = 0; i < xRes * yRes; ++i)
        hrgba[i] = Rgba(pixels[3*i], pixels[3*i+1], pixels[3*i+2],
                        alpha ? alpha[i]: 1.f);

    Box2i displayWindow(V2i(0,0), V2i(totalXRes-1, totalYRes-1));
    Box2i dataWindow(V2i(xOffset, yOffset), V2i(xOffset + xRes - 1, yOffset + yRes - 1));

    try {
        RgbaOutputFile file(name.c_str(), displayWindow, dataWindow, WRITE_RGBA);
        file.setFrameBuffer(hrgba - xOffset - yOffset * xRes, 1, xRes);
        file.writePixels(yRes);
    }
    catch (const std::exception &e) {
        Error("Unable to write image file \"%s\": %s", name.c_str(),
            e.what());
    }
}
