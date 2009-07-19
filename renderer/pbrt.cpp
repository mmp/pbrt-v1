
/*
 * pbrt source code Copyright(c) 1998-2007 Matt Pharr and Greg Humphreys
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 * (See file License.txt for complete license)
 */

// pbrt.cpp*
#include "pbrt.h"
#include "api.h"
// main program
int main(int argc, char *argv[]) {
	// Print welcome banner
	printf("pbrt version %1.3f of %s at %s\n",
	       PBRT_VERSION, __DATE__, __TIME__);
	printf("Copyright (c)1998-2007 Matt Pharr and "
	       "Greg Humphreys.\n");
	printf("The source code to pbrt (but *not* the contents of the book) is\n");
	printf("covered by the GNU General Public License.  See the file COPYING.txt\n");
	printf("for the conditions of the license.\n");
	fflush(stdout);
	pbrtInit();
	// Process scene description
	if (argc == 1) {
		// Parse scene from standard input
		ParseFile("-");
	} else {
		// Parse scene from input files
		for (int i = 1; i < argc; i++)
			if (!ParseFile(argv[i]))
				Error("Couldn't open scene file \"%s\"\n", argv[i]);
	}
	pbrtCleanup();
	return 0;
}
