
CC=			gcc
CXX=		g++
CFLAGS=		-Wall -g -O2  
CXXFLAGS=	$(CFLAGS)
DFLAGS=		-DFILE_OFFSET_BITS=64#-DHAVE_PTHREAD #-DDEBUG#-D_FILE_OFFSET_BITS=64 
SRC= 		utils.c malloc_wrap.c bwt.c bwtio.c is.c\
			kstring.c bntseq.c bwtmisc.c\
		bwt_gen.c QSufSort.c\
		bitmap.c variants.c lookup.c query.c setFileName.c
PROG=		index aln gen_seed aln_debug test_index test_aln
INCLUDES=	
LIBS=		-lm -lz
SUBDIRS=	. 

all: index 
gen_seed:
		$(CC) $(CXXFLAGS) $(DFLAGS) $(SRC) index1.c gen_test_seeds.c -o $@ $(LIBS)
		mv gen_seed Bin/
index: 
		$(CC) $(CFLAGS) $(DFLAGS) -DMAIN_INDEX $(SRC)  index1.c -o $@ $(LIBS)
		mv index ../Bin/
aln:
		$(CC) $(CFLAGS) $(DFLAGS) -DALN_MAIN $(SRC) index1.c bsearch.c aln.c -o $@ $(LIBS)
		mv aln Bin/
test_aln:
		$(CC) $(CFLAGS) $(DFLAGS) -DALN_MAIN $(SRC) index1.c bsearch.c test_aln.c -o $@ $(LIBS)
		mv test_aln Bin/
aln_debug:
		$(CC) $(CFLAGS) $(DFLAGS) -DALN_MAIN $(SRC) index1.c bsearch.c aln.c -o $@ $(LIBS)
		mv aln_debug Debug/
test_index:
		$(CC) $(CFLAGS) $(DFLAGS) -DMAIN_INDEX $(SRC) test_index.c -o $@ $(LIBS)
		mv test_index Bin/

clean:
		rm -f gmon.out *.o a.out $(PROG) *~ *.a
		rm -f Bin/*

