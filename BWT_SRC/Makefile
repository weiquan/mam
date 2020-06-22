
CC=			gcc
CXX=		g++
CFLAGS=		-Wall -g  
CXXFLAGS=	$(CFLAGS)
DFLAGS=		-DFILE_OFFSET_BITS=64#-DHAVE_PTHREAD #-DDEBUG#-D_FILE_OFFSET_BITS=64 
SRC= 		utils.c malloc_wrap.c bwt.c bwtio.c is.c\
			kstring.c bntseq.c bwtmisc.c\
		bwt_gen.c QSufSort.c\
		bitmap.c lookup.c query.c setFileName.c rbwt_alt.c
PROG=		index aln gen_seed aln_debug test_index test_aln
INCLUDES=	
LIBS=		-lm -lz
SUBDIRS=	. 

all: index 
index: 
		$(CC) $(CFLAGS) -O2 $(DFLAGS) -DMAIN_INDEX $(SRC)  build_hier_idx.c index.c -o $@ $(LIBS)
		mv $@ ../Bin/
aln: 
		$(CC) $(CFLAGS) -O2 $(DFLAGS) -DMAIN_ALN $(SRC)  index.c build_hier_idx.c aln.c -o $@ $(LIBS)
		mv $@ ../Bin/

index_debug: 
		$(CC) -Wall -g $(DFLAGS) -DMAIN_INDEX $(SRC)  build_hier_idx.c index.c -o $@ $(LIBS)
		mv $@ ../Bin/
aln_debug: 
		$(CC) -Wall -g $(CFLAGS) $(DFLAGS) -DMAIN_ALN -DDEBUG  $(SRC) index.c build_hier_idx.c aln.c -o $@ $(LIBS)
		mv $@ ../Bin/



libabwt: 
		$(CC) -c $(DFLAGS) $(SRC)  index.c ${LIBS}
		ar crv libabwt.a *.o
		rm *.o
		#mv libabwt.a test
clean:
		rm -f gmon.out *.o a.out $(PROG) *~ *.a
		rm -f Bin/*

