
CC=			gcc
CXX=		g++
CFLAGS=		-Wall -g #-O2  
CXXFLAGS=	$(CFLAGS)
DFLAGS=		-DFILE_OFFSET_BITS=64#-DHAVE_PTHREAD #-DDEBUG#-D_FILE_OFFSET_BITS=64 
SRC= 		SeedExtAlign2.c lookup.c ksw.c query.c SeedExtAlign.c	



PROG=		aln aln_debug 
INCLUDES=	
LIBS=		-lm -lz
SUBDIRS=	. 

all:$(PROG)
aln:
		$(CC) $(CFLAGS) -O2 $(DFLAGS) $(SRC)  -o $@ $(LIBS)
		mv $@ ../../Index/
aln_debug:
		$(CC) $(CFLAGS) $(DFLAGS) $(SRC)  -o $@ $(LIBS)


