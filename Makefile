# user may define these whole flags
# LDFLAGS
# CPPFLAGS
# CXXFLAGS
# CFLAGS

# or these user-set flags that will be added to standard flags
CXXFLAGSextra ?=

# user may define the compiler
CXX ?= 

# pre-defined flags

COMPTIMEPLACE := -D'COMPILATION_TIME_PLACE="/dev/null"'

CXXFLAGS_common := -pipe -std=c++11 -Wall -Wextra -fPIC $(COMPTIMEPLACE)
CXXFLAGS_main := -O3 -g $(CXXFLAGS_common)
CXXFLAGS_gdb :=  -O0 -g $(CXXFLAGS_common)

CFLAGS := -O3 -pipe -Wall -Wextra -fPIC $(CFLAGS)


##########################################################################################################

OBJECTS = ParametersChimeric_initialize.o ParametersSolo.o SoloRead.o SoloRead_record.o \
	SoloReadBarcode.o SoloReadBarcode_getCBandUMI.o \
	SoloReadFeature.o SoloReadFeature_record.o SoloReadFeature_inputRecords.o \
	Solo.o SoloFeature.o SoloFeature_collapseUMI.o SoloFeature_outputResults.o SoloFeature_processRecords.o\
	ReadAlign_outputAlignments.o  \
	ReadAlign.o \
	SharedMemory.o PackedArray.o SuffixArrayFuns.o Parameters.o InOutStreams.o SequenceFuns.o Genome.o Stats.o \
	Transcript.o Transcript_alignScore.o Transcript_generateCigarP.o Chain.o \
	Transcript_variationAdjust.o Variation.o ReadAlign_waspMap.o \
	ReadAlign_storeAligns.o ReadAlign_stitchPieces.o ReadAlign_multMapSelect.o ReadAlign_mapOneRead.o readLoad.o \
	ReadAlignChunk.o \
	OutSJ.o outputSJ.o blocksOverlap.o ThreadControl.o sysRemoveDir.o \
	ReadAlign_maxMappableLength2strands.o binarySearch2.o\
	ReadAlign_outputTranscriptSAM.o ReadAlign_outputTranscriptSJ.o ReadAlign_outputTranscriptCIGARp.o \
	ReadAlign_createExtendWindowsWithAlign.o ReadAlign_assignAlignToWindow.o ReadAlign_oneRead.o \
	ReadAlign_stitchWindowSeeds.o \
	ReadAlign_peOverlapMergeMap.o ReadAlign_mappedFilter.o \
	ReadAlign_chimericDetection.o ReadAlign_chimericDetectionOld.o ReadAlign_chimericDetectionOldOutput.o\
	ChimericDetection.o ChimericDetection_chimericDetectionMult.o ReadAlign_chimericDetectionPEmerged.o \
	ChimericAlign.o ChimericSegment.o ReadAlign_calcCIGAR.o \
	ChimericAlign_chimericJunctionOutput.o ChimericAlign_chimericStitching.o \
	stitchWindowAligns.o extendAlign.o stitchAlignToTranscript.o alignSmithWaterman.o \
	Genome_genomeGenerate.o genomeParametersWrite.o genomeScanFastaFiles.o genomeSAindex.o \
	Genome_insertSequences.o insertSeqSA.o funCompareUintAndSuffixes.o funCompareUintAndSuffixesMemcmp.o \
	TimeFunctions.o ErrorWarning.o loadGTF.o streamFuns.o stringSubstituteAll.o \
	Transcriptome.o Transcriptome_quantAlign.o Transcriptome_geneFullAlignOverlap.o \
	ReadAlign_quantTranscriptome.o Quantifications.o Transcriptome_geneCountsAddAlign.o \
	sjdbLoadFromFiles.o sjdbLoadFromStream.o sjdbPrepare.o sjdbBuildIndex.o sjdbInsertJunctions.o \
	Parameters_openReadsFiles.o Parameters_closeReadsFiles.o Parameters_readSAMheader.o \
	bam_cat.o serviceFuns.o GlobalVariables.o \
	BAMoutput.o BAMfunctions.o ReadAlign_alignBAM.o BAMbinSortByCoordinate.o BAMbinSortUnmapped.o \
	orbit.o

SOURCES := $(wildcard *.cpp) $(wildcard *.c)


%.o : %.cpp
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $<

%.o : %.c
	$(CXX) -c $(CPPFLAGS) $(CFLAGS) $<

all: STAR

.PHONY: clean
clean:
	rm -f *.o *.a Depend.list

.PHONY: CLEAN
CLEAN:
	rm -f *.o STAR Depend.list
	$(MAKE) -C htslib clean

.PHONY: cleanRelease
cleanRelease:
	rm -f *.o Depend.list
	$(MAKE) -C htslib clean

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),cleanRelease)
ifneq ($(MAKECMDGOALS),CLEAN)
Depend.list: $(SOURCES) parametersDefault.xxd htslib
	echo $(SOURCES)
	'rm' -f ./Depend.list
	$(CXX) $(CXXFLAGS_common) -MM $^ >> Depend.list
include Depend.list
endif
endif
endif

htslib : htslib/libhts.a

htslib/libhts.a :
	$(MAKE) -C htslib lib-static

parametersDefault.xxd: parametersDefault
	xxd -i parametersDefault > parametersDefault.xxd

liborbit.a : CXXFLAGS := $(CXXFLAGSextra) $(CXXFLAGS_main) $(CXXFLAGS)
liborbit.a : $(OBJECTS)
	ar -csru $@ $(OBJECTS)


