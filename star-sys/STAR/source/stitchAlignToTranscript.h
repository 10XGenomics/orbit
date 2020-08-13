#include "IncludeDefine.h"
#include "Parameters.h"
#include "Transcript.h"
#include "Genome.h"

intScore stitchAlignToTranscript(uint rAend, uint gAend, uint rBstart, uint gBstart, uint L, uint iFragB, uint sjAB, const Parameters& P, char* R, const Genome &mapGen,  Transcript *trA, uint outFilterMismatchNmaxTotal);

