#include "IncludeDefine.h"
#include "Parameters.h"
#include "Transcript.h"
#include "ReadAlign.h"
#include "BAMfunctions.h"
#include "blocksOverlap.h"

void ReadAlign::chimericDetection() {

    chimRecord=false;

    if (P.pCh.segmentMin==0) {//no chimeric detection requested
        return;
    };
    if (P.outFilterBySJoutStage>1) {//no chimeric output for stage=2. REVISIT: NOT SURE why
        return;
    };

    //output chains for out-of-STAR chimeric detection
    #ifdef OUTPUT_localChains
    {
        P.inOut->outLocalChains << readName <<"\t"<< Read0[0] <<"\t"<< Read0[1] << "\n";
        for (uint iw=0; iw<nW; iw++) {
            for (uint itr=0;itr<nWinTr[iw];itr++) {
                P.inOut->outLocalChains << trAll[iw][itr]->maxScore<<"\t"<< trAll[iw][itr]->Chr<<"\t"<<trAll[iw][itr]->Str<<"\t"<<trAll[iw][itr]->nExons;
                for (uint ib=0;ib<trAll[iw][itr]->nExons;ib++) {
                    P.inOut->outLocalChains <<"\t"<< trAll[iw][itr]->exons[ib].G-mapGen.chrStart[trAll[iw][itr]->Chr] \
                                             <<"\t"<< trAll[iw][itr]->exons[ib].R <<"\t"<< trAll[iw][itr]->exons[ib].L;
                };
                P.inOut->outLocalChains <<"\n";
            };
        };
    };
    #endif


    if (P.pCh.multimapNmax==0) {
        chimRecord=chimericDetectionOld();
        chimericDetectionOldOutput();
    } else if (trBest->maxScore <= (int) (readLength[0]+readLength[1]) - (int) P.pCh.nonchimScoreDropMin) {//require big enough drop in the best score
        chimRecord=chimDet->chimericDetectionMult(nW, readLength, trBest->maxScore, false);
    };

    if ( chimRecord ) {
        statsRA.chimericAll++;
    };

    return;
};//END
