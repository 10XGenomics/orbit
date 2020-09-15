#include "ChimericDetection.h"

ChimericDetection::ChimericDetection(const Parameters &Pin, Transcript ***trAll, uint *nWinTr, char** Read1in, const Genome &mapGenIn, fstream *ostreamChimJunctionIn, ReadAlign *RAin)
                  : P(Pin), RA(RAin), trAll(trAll), nWinTr(nWinTr), Read1(Read1in), outGen(mapGenIn), ostreamChimJunction(ostreamChimJunctionIn)
                                    {
};
