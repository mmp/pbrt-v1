#ifndef BISON_CORE_PBRTPARSE_HPP
# define BISON_CORE_PBRTPARSE_HPP

#ifndef YYSTYPE
typedef union {
char string[1024];
float num;
ParamArray *ribarray;
} yystype;
# define YYSTYPE yystype
# define YYSTYPE_IS_TRIVIAL 1
#endif
# define	STRING	257
# define	ID	258
# define	NUM	259
# define	LBRACK	260
# define	RBRACK	261
# define	ACCELERATOR	262
# define	AREALIGHTSOURCE	263
# define	ATTRIBUTEBEGIN	264
# define	ATTRIBUTEEND	265
# define	CAMERA	266
# define	CONCATTRANSFORM	267
# define	COORDINATESYSTEM	268
# define	COORDSYSTRANSFORM	269
# define	FILM	270
# define	IDENTITY	271
# define	LIGHTSOURCE	272
# define	LOOKAT	273
# define	MATERIAL	274
# define	OBJECTBEGIN	275
# define	OBJECTEND	276
# define	OBJECTINSTANCE	277
# define	PIXELFILTER	278
# define	REVERSEORIENTATION	279
# define	ROTATE	280
# define	SAMPLER	281
# define	SCALE	282
# define	SEARCHPATH	283
# define	SHAPE	284
# define	SURFACEINTEGRATOR	285
# define	TEXTURE	286
# define	TRANSFORMBEGIN	287
# define	TRANSFORMEND	288
# define	TRANSFORM	289
# define	TRANSLATE	290
# define	VOLUME	291
# define	VOLUMEINTEGRATOR	292
# define	WORLDBEGIN	293
# define	WORLDEND	294
# define	HIGH_PRECEDENCE	295


extern YYSTYPE yylval;

#endif /* not BISON_CORE_PBRTPARSE_HPP */
