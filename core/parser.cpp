
/*
 * pbrt source code Copyright(c) 1998-2007 Matt Pharr and Greg Humphreys
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 * (See file License.txt for complete license)
 */

// parser.cpp*
#include "pbrt.h"
// Parsing Global Interface
COREDLL bool ParseFile(const char *filename) {
	extern FILE *yyin;
	extern int yyparse(void);
	extern string current_file;
	extern int line_num;
	extern int yydebug;

	if (getenv("PBRT_YYDEBUG") != NULL)
		yydebug = 1;

	if (strcmp(filename, "-") == 0)
		yyin = stdin;
	else
		yyin = fopen(filename, "r");
	if (yyin != NULL) {
		current_file = filename;
		if (yyin == stdin) current_file = "<standard input>";
		line_num = 1;
		yyparse();
		if (yyin != stdin) fclose(yyin);
	}
	current_file = "";
	line_num = 0;
	return (yyin != NULL);
}
