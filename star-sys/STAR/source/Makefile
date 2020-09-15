# user may define these whole flags
# CPPFLAGS
# CXXFLAGS
UNAME := $(shell uname)

# or these user-set flags that will be added to standard flags
CXXFLAGSextra ?=

# user may define the compiler
CXX ?= 

# pre-defined flags

COMPTIMEPLACE := -D'COMPILATION_TIME_PLACE="source/Makefile"'

ifeq ($(UNAME), Darwin)
	CXXFLAGS_extra := -D'COMPILE_FOR_MAC'
else
	CXXFLAGS_extra :=
endif
CXXFLAGS_common := -pipe -std=c++11 -Wall -Wextra -Werror -fPIC $(CXXFLAGS_extra) $(COMPTIMEPLACE)

CXXFLAGS_main := -O3 -g $(CXXFLAGS_common)
CXXFLAGS_gdb :=  -O0 -g $(CXXFLAGS_common)


##########################################################################################################

OBJECTS = orbit.o \
    InOutStreams.o \
    Parameters.o \
    ParametersSolo.o \
    ParametersChimeric_initialize.o \
    PackedArray.o \
    SequenceFuns.o \
    Genome.o \
    Genome_insertSequences.o \
    Genome_genomeGenerate.o \
    streamFuns.o \
    genomeScanFastaFiles.o \
    TimeFunctions.o \
    insertSeqSA.o \
    ReadAlign.o \
    Transcript.o \
    Transcriptome_quantAlign.o \
    funCompareUintAndSuffixesMemcmp.o \
    genomeSAindex.o \
    ReadAlign_outputAlignments.o \
    ReadAlign_outputTranscriptSAM.o \
    ReadAlign_quantTranscriptome.o \
    ReadAlign_calcCIGAR.o \
    ReadAlign_storeAligns.o \
    SuffixArrayFuns.o \
    ReadAlign_oneRead.o \
    ReadAlign_mapOneRead.o \
    ReadAlign_stitchPieces.o \
    ReadAlign_mappedFilter.o \
    ReadAlign_maxMappableLength2strands.o \
    ReadAlign_assignAlignToWindow.o \
    ReadAlign_createExtendWindowsWithAlign.o \
    ReadAlign_multMapSelect.o \
    readLoad.o \
    stitchWindowAligns.o \
    extendAlign.o \
    binarySearch2.o \
    blocksOverlap.o \
    stitchAlignToTranscript.o \
    ErrorWarning.o

%.o : %.cpp
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $<

all: liborbit.a

clean:
	rm -f *.o *.a

parametersDefault.xxd: parametersDefault
	xxd -i parametersDefault > parametersDefault.xxd

liborbit.a : CXXFLAGS := $(CXXFLAGSextra) $(CXXFLAGS_main) $(CXXFLAGS)
liborbit.a : $(OBJECTS)
	ar -csru $@ $(OBJECTS)


