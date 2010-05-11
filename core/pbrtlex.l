
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

/* state used for include file stuff */
%{

#define YY_MAIN 0
#define YY_NEVER_INTERACTIVE 1

#include "pbrt.h"
#include "api.h"

struct ParamArray;

#ifdef WIN32
#include "pbrtparse.hpp"
#else
#include "pbrtparse.h"
#endif

#ifdef WIN32
#pragma warning ( disable: 4244 )
#endif

struct IncludeInfo {
	string filename;
	YY_BUFFER_STATE bufState;
	int lineNum;
};
vector<IncludeInfo> includeStack;

extern int line_num;
int str_pos;

void add_string_char( char c )
{
	yylval.string[str_pos++] = c;
	yylval.string[str_pos] = '\0';
}

extern void yyerror( const char *str );
extern void ParseError( const char *format, ... ) PRINTF_FUNC;

void include_push(char *filename) {
	if (includeStack.size() > 32)
		Severe("Only 32 levels of nested Include allowed in scene files.");
	IncludeInfo ii;
	extern string current_file;
	ii.filename = current_file;
	ii.bufState = YY_CURRENT_BUFFER;
	ii.lineNum = line_num;
	includeStack.push_back(ii);

	current_file = filename;
	line_num = 1;

	yyin = fopen(filename, "r");
	if (!yyin)
		Severe("Unable to open included scene file \"%s\"", filename);
	yy_switch_to_buffer(yy_create_buffer(yyin, YY_BUF_SIZE));
}

void include_pop(void)
{
	extern int line_num;
	extern string current_file;
	yy_delete_buffer(YY_CURRENT_BUFFER);
	yy_switch_to_buffer(includeStack.back().bufState);
	current_file = includeStack.back().filename;
	line_num = includeStack.back().lineNum;
	includeStack.pop_back();
}
%}
%option nounput
WHITESPACE [ \t\r]+
NUMBER [-+]?([0-9]+|(([0-9]+\.[0-9]*)|(\.[0-9]+)))([eE][-+]?[0-9]+)?
IDENT [a-zA-Z_][a-zA-Z_0-9]*
%x STR COMMENT INCL INCL_FILE
%%
Include { BEGIN(INCL); }
<INCL>{WHITESPACE} ;
<INCL>\" { BEGIN(INCL_FILE); }
<INCL>. { ParseError( "Illegal character following Include directive" ); }
<INCL_FILE>\" { BEGIN INITIAL; }
<INCL_FILE>. { ParseError( "Illegal character in Include file name" ); }
<INCL_FILE>[\-a-zA-Z_\.0-9/ \t]+ {
	BEGIN(INITIAL);
	include_push(yytext);
}
"#" { BEGIN COMMENT; }
<COMMENT>. /* eat it up */
<COMMENT>\n { line_num++; BEGIN INITIAL; }
Accelerator		{ return ACCELERATOR; }
AreaLightSource		{ return AREALIGHTSOURCE;	}
AttributeBegin 		{ return ATTRIBUTEBEGIN; 	}
AttributeEnd 		{ return ATTRIBUTEEND; 		}
Camera			{ return CAMERA;		}
ConcatTransform 	{ return CONCATTRANSFORM; 	}
CoordinateSystem	{ return COORDINATESYSTEM;	}
CoordSysTransform	{ return COORDSYSTRANSFORM;	}
Film	 		{ return FILM;	 		}
Identity		{ return IDENTITY;		}
LightSource 		{ return LIGHTSOURCE;		}
LookAt			{ return LOOKAT;		}
Material		{ return MATERIAL;		}
ObjectBegin 		{ return OBJECTBEGIN;		}
ObjectEnd 		{ return OBJECTEND;		}
ObjectInstance 		{ return OBJECTINSTANCE;	}
PixelFilter 		{ return PIXELFILTER; 		}
ReverseOrientation 	{ return REVERSEORIENTATION;	}
Rotate			{ return ROTATE;		}
Sampler			{ return SAMPLER;		}
SearchPath		{ return SEARCHPATH;		}
Scale 			{ return SCALE;			}
Shape			{ return SHAPE;			}
SurfaceIntegrator	{ return SURFACEINTEGRATOR;	}
Texture			{ return TEXTURE;		}
TransformBegin	{ return TRANSFORMBEGIN;	}
TransformEnd	{ return TRANSFORMEND;	}
Transform		{ return TRANSFORM;		}
Translate		{ return TRANSLATE; 		}
Volume			{ return VOLUME;		}
VolumeIntegrator	{ return VOLUMEINTEGRATOR;	}
WorldBegin 		{ return WORLDBEGIN; 		}
WorldEnd 		{ return WORLDEND;		}
{WHITESPACE} /* do nothing */
\n { line_num++; }
{NUMBER} {
  yylval.num = (float) atof(yytext);
  return NUM;
}
{IDENT} {
	strcpy( yylval.string, yytext );
	return ID;
}
"[" { return LBRACK; }
"]" { return RBRACK; }
\" { BEGIN STR; str_pos = 0; }
<STR>\\n {add_string_char('\n');}
<STR>\\t {add_string_char('\t');}
<STR>\\r {add_string_char('\r');}
<STR>\\b {add_string_char('\b');}
<STR>\\f {add_string_char('\f');}
<STR>\\\" {add_string_char('\"');}
<STR>\\\\ {add_string_char('\\');}
<STR>\\[0-9]{3} {
  int val = atoi(yytext+1);
  while(val > 256)
    val -= 256;
  add_string_char(val);
}
<STR>\\\n {line_num++;}
<STR>\\. { add_string_char(yytext[1]);}
<STR>\" {BEGIN INITIAL; return STRING;}
<STR>. {add_string_char(yytext[0]);}
<STR>\n {ParseError("Unterminated string!");}

. { ParseError( "Illegal character: %c",yytext[0] ); }
%%
int yywrap(void)
{
	if (includeStack.size() ==0) return 1;
	include_pop();
	BEGIN(INCL_FILE);
	return 0;
}

