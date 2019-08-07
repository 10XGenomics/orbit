#ifndef CODE_ChimericDetection
#define CODE_ChimericDetection

#include "IncludeDefine.h"
#include "Parameters.h"
#include "Transcript.h"
#include "ChimericAlign.h"
#include "Genome.h"

class ReadAlign;

class ChimericDetection {
    private:
        const Parameters &P;
        ReadAlign *RA;
        Transcript ***trAll;
        uint nW, *nWinTr;
        char** Read1;
        const Genome &outGen;
        uint *readLength;

    public:
        uint chimN;
        vector <ChimericAlign> chimAligns;
        bool chimRecord;
        int chimScoreBest;

        ChimericDetection(const Parameters &Pin, Transcript ***trAll, uint *nWinTr, char** Read1in, const Genome &genomeIn, fstream *ostreamChimJunctionIn, ReadAlign *RA);
        bool chimericDetectionMult(uint nWin, uint *readLengthIn, int maxNonChimAlignScore, bool PEmerged_flag);
        bool chimericDetectionMult(uint nWin, uint *readLengthIn);
        fstream *ostreamChimJunction;
};

#endif
