/*
 g++ -g -Wall exrassemble.cpp -I../../src/tangled_code/OpenEXR/include -ltiff -L../../src/tangled_code/OpenEXR/lib-linux -lIlmImf -lImath -lIex -lHalf  

*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/types.h>
#include <dirent.h>
#include <string.h>
#include <ImfAttribute.h>
#include <ImfInputFile.h>
#include <ImfOutputFile.h>
#include <ImfChannelList.h>
#include <ImfFrameBuffer.h>
#include <half.h>

using namespace Imf;
using namespace Imath;

void WriteEXR(const char *name, half *rgba, int xRes, int yRes);
half *ReadEXR(const char *name, int *x0, int *y0, int *xSize, int *ySize,
	      int *xres, int *yres);

int main(int argc, char *argv[]) 
{
    if (argc != 3) {
	fprintf(stderr, "usage: exrtogether [exr dir] [foo.exr]\n");
	return 1;
    }

    const char *dirname = argv[1];
    const char *outname = argv[2];

    int xres, yres;
    int ndone = 0;
    half *rgba = 0;
	
    DIR *dir = opendir(dirname);
    struct dirent *d;
    while ((d = readdir(dir))) {
	if (!strstr(d->d_name, ".exr")) continue;
	int x0, y0, xr, yr, totx, toty;
	char fn[1024];
	sprintf(fn, "%s/%s", dirname, d->d_name);
	half *subpix = ReadEXR(fn, &x0, &y0, &xr, &yr, &totx, &toty);
	if (!subpix) {
	    fprintf(stderr, "couldn't read exr file \"%s\"!\n", fn);
	    continue;
	}

	if (!rgba) {
	    xres = totx;
	    yres = toty;
	    rgba = new half[xres*yres*4];
	    memset(rgba, 0, xres*yres*4*sizeof(half));
	}
	else
	    assert(xres == totx && yres == toty);

	ndone += xr*yr;
	for (int i = 0; i < yr; ++i)
	    memcpy(rgba + 4 * (x0 + xres * (y0 + i)),
		   subpix + 4 * i * xr, 4 * xr * sizeof(half));
	delete[] subpix;
    }
    closedir(dir);

    fprintf(stderr, "Got %f%% of image\n", 100.f * ndone / (xres *yres));
    WriteEXR(outname, rgba, xres, yres);

    return 0;
}

void WriteEXR(const char *name, half *rgba, int xRes, int yRes) 
{
    Header header(xRes, yRes);
    header.channels().insert("R", Channel (HALF));
    header.channels().insert("G", Channel (HALF));
    header.channels().insert("B", Channel (HALF));
    header.channels().insert("A", Channel (HALF));

    int stride = 4;

    FrameBuffer fb;
    fb.insert("R", Slice(HALF, (char *)rgba, stride*sizeof(half),
			 stride*xRes*sizeof(half)));
    fb.insert("G", Slice(HALF, (char *)rgba+sizeof(half), stride*sizeof(half),
			 stride*xRes*sizeof(half)));
    fb.insert("B", Slice(HALF, (char *)rgba+2*sizeof(half), stride*sizeof(half),
			 stride*xRes*sizeof(half)));
    fb.insert("A", Slice(HALF, (char *)rgba+3*sizeof(half), stride*sizeof(half),
			 stride*xRes*sizeof(half)));

    OutputFile file(name, header);
    file.setFrameBuffer(fb);
    file.writePixels(yRes);
}

half *ReadEXR(const char *name, int *x0, int *y0, int *xSize, int *ySize,
	      int *totxres, int *totyres)
{
    InputFile file(name);
    Box2i dw = file.header().dataWindow();
    Box2i totalImage = file.header().displayWindow();

    *totxres = totalImage.max.x - totalImage.min.x + 1;
    *totyres = totalImage.max.y - totalImage.min.y + 1;
    *x0 = dw.min.x;
    *y0 = dw.min.y;
    *xSize = dw.max.x - dw.min.x + 1;
    *ySize = dw.max.y - dw.min.y + 1;

    half *rgba = new half[4 * *xSize * *ySize];

    half *rstart = rgba - (4 * (*x0 + *y0 * *xSize));

    FrameBuffer frameBuffer;
    frameBuffer.insert("R", Slice(HALF, (char *)rstart,
				  4*sizeof(half), *xSize * 4 * sizeof(half), 1, 1, 0.0));
    frameBuffer.insert("G", Slice(HALF, (char *)rstart+sizeof(half),
				  4*sizeof(half), *xSize * 4 * sizeof(half), 1, 1, 0.0));
    frameBuffer.insert("B", Slice(HALF, (char *)rstart+2*sizeof(half),
				  4*sizeof(half), *xSize * 4 * sizeof(half), 1, 1, 0.0));
    frameBuffer.insert("A", Slice(HALF, (char *)rstart+3*sizeof(half),
				  4*sizeof(half), *xSize * 4 * sizeof(half), 1, 1, 1.0));

    file.setFrameBuffer(frameBuffer);
    file.readPixels(dw.min.y, dw.max.y);

    return rgba;
}
