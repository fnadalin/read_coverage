
CCC := gcc
CPP := g++
CPPFLAGS := -I samtools -O3 -std=c++11
CFLAGS := -I samtools -O3 -std=c++0x -fpermissive
LDFLAGS := -lz

DEPS := bam.h bam_aux.h bam_import.h bam_index.h bam_pileup.h bgzf.h faidx.h knetfile.h \
 kstring.h razf.h sam.h sam_header.h kaln.h klist.h kstring.h bam_endian.h common.h

OBJ := Sam.o Bed.o samtools/bam.o samtools/bam_aux.o samtools/bam_import.o samtools/bam_index.o samtools/bam_pileup.o samtools/bgzf.o samtools/faidx.o \
 samtools/knetfile.o samtools/kstring.o samtools/razf.o samtools/sam.o samtools/sam_header.o read-coverage.o

OBJ1 := Sam.o Bed.o samtools/bam.o samtools/bam_aux.o samtools/bam_import.o samtools/bam_index.o samtools/bam_pileup.o samtools/bgzf.o samtools/faidx.o \
 samtools/knetfile.o samtools/kstring.o samtools/razf.o samtools/sam.o samtools/sam_header.o sj-coverage.o

OBJ2 := Sam.o Bed.o samtools/bam.o samtools/bam_aux.o samtools/bam_import.o samtools/bam_index.o samtools/bam_pileup.o samtools/bgzf.o samtools/faidx.o \
 samtools/knetfile.o samtools/kstring.o samtools/razf.o samtools/sam.o samtools/sam_header.o repeat-coverage-single-cell.o

samtools/%.o: samtools/%.c $(DEPS)
	$(CCC) $(CFLAGS) -c $< -o $@

%.o: %.cpp $(DEPS)
	$(CPP) $(CPPFLAGS) -c $< -o $@

all: read-coverage sj-coverage repeat-coverage-single-cell

read-coverage: $(OBJ)
	$(CPP) $(CPPFLAGS) $(OBJ) -o $@ $(LDFLAGS)

sj-coverage: $(OBJ1)
	$(CPP) $(CPPFLAGS) $(OBJ1) -o $@ $(LDFLAGS)

repeat-coverage-single-cell: $(OBJ2)
	$(CPP) $(CPPFLAGS) $(OBJ2) -o $@ $(LDFLAGS)


