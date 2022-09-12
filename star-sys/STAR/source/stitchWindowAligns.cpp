#include "stitchWindowAligns.h"
#include "blocksOverlap.h"
#include "ErrorWarning.h"
#include "binarySearch2.h"
#include <cmath>
#include <vector>

void stitchWindowAligns(uint iA, int Score, bool WAincl[], uint tR2, uint tG2, Transcript trA, \
                        uint Lread, vector<uiWA>& WA, char* R, const Genome &mapGen, \
                        const Parameters& P, Transcript** wTr, uint* nWinTr, ReadAlign *RA) {
    //recursively stitch aligns for one gene
    //*nWinTr - number of transcripts for the current window
    uint nA = WA.size();
    if (iA>=nA && tR2==0) return; //no aligns in the transcript

    if (iA>=nA) {//no more aligns to add, finalize the transcript

        //extend first
        Transcript trAstep1;

        int vOrder[2]; //decide in which order to extend: extend the 5' of the read first

        #if EXTEND_ORDER==1
        if ( trA.roStr==0 ) {//decide in which order to extend: extend the 5' of the read first
            vOrder[0]=0; vOrder[1]=1;
        } else {
            vOrder[0]=1; vOrder[1]=0;
        };
        #elif EXTEND_ORDER==2
            vOrder[0]=0; vOrder[1]=1;
        #else
            #error "EXTEND_ORDER value unrecognized"
        #endif

        for (int iOrd=0;iOrd<2;iOrd++) {

            switch (vOrder[iOrd]) {

            case 0: //extend at start

            if (trA.rStart>0) {// if transcript does not start at base, extend to the read start
                trAstep1.reset();
                uint imate=trA.exons[0].iFrag;
                if ( extendAlign(R, mapGen.G, trA.rStart-1, trA.gStart-1, -1, -1, trA.rStart, tR2-trA.rStart+1, \
                                 trA.nMM, RA->outFilterMismatchNmaxTotal, P.outFilterMismatchNoverLmax, \
                                 P.alignEndsType.ext[imate][(int)(trA.Str!=imate)], &trAstep1) ) {//if could extend

                    trA.add(&trAstep1);
                    Score += trAstep1.maxScore;

                    trA.exons[0].R = trA.rStart = trA.rStart - trAstep1.extendL;
                    trA.exons[0].G = trA.gStart = trA.gStart - trAstep1.extendL;
                    trA.exons[0].L += trAstep1.extendL;

                };
            //TODO penalize the unmapped bases at the start
            };
            break;

            case 1: //extend at end

            if ( tR2<Lread ) {//extend alignment to the read end
                trAstep1.reset();
                uint imate=trA.exons[trA.nExons-1].iFrag;
                if ( extendAlign(R, mapGen.G, tR2+1, tG2+1, +1, +1, Lread-tR2-1, tR2-trA.rStart+1, \
                                 trA.nMM, RA->outFilterMismatchNmaxTotal,  P.outFilterMismatchNoverLmax, \
                                 P.alignEndsType.ext[imate][(int)(imate==trA.Str)], &trAstep1) ) {//if could extend

                    trA.add(&trAstep1);
                    Score += trAstep1.maxScore;

                    tR2 += trAstep1.extendL;
                    tG2 += trAstep1.extendL;

                    trA.exons[trA.nExons-1].L += trAstep1.extendL;//extend the length of the last exon

                };
            //TODO penalize unmapped bases at the end
            };
        };
        };

        if (!P.alignSoftClipAtReferenceEnds.yes &&  \
                ( (trA.exons[trA.nExons-1].G + Lread-trA.exons[trA.nExons-1].R) > (mapGen.chrStart[trA.Chr]+mapGen.chrLength[trA.Chr]) || \
                   trA.exons[0].G<(mapGen.chrStart[trA.Chr]+trA.exons[0].R) ) ) {
            return; //no soft clipping past the ends of the chromosome
        };


        trA.rLength = 0;
        for (uint isj=0;isj<trA.nExons;isj++) {
            trA.rLength += trA.exons[isj].L;
        };
        trA.gLength = tG2+1-trA.gStart;

        //check exons lengths including repeats, do not report a transcript with short exons
        for (uint isj=0;isj<trA.nExons-1;isj++) {//check exons for min length, if they are not annotated and precede a junction
            if ( trA.canonSJ[isj]>=0 ) {//junction
                if (trA.sjAnnot[isj]==1) {//sjdb
                    if (  ( trA.exons[isj].L   < P.alignSJDBoverhangMin && (isj==0            || trA.canonSJ[isj-1]==-3 || (trA.sjAnnot[isj-1]==0 && trA.canonSJ[isj-1]>=0) ) )\
                       || ( trA.exons[isj+1].L < P.alignSJDBoverhangMin && (isj==trA.nExons-2 || trA.canonSJ[isj+1]==-3 || (trA.sjAnnot[isj+1]==0 && trA.canonSJ[isj+1]>=0) ) ) )return;
                } else {//non-sjdb
                    if (  trA.exons[isj].L < P.alignSJoverhangMin + trA.shiftSJ[isj][0] \
                       || trA.exons[isj+1].L < P.alignSJoverhangMin + trA.shiftSJ[isj][1]   ) return;
                };
            };
        };
        if (trA.nExons>1 && trA.sjAnnot[trA.nExons-2]==1 && trA.exons[trA.nExons-1].L < P.alignSJDBoverhangMin) return; //this exon was not checkedin the cycle above

        //filter strand consistency
        uint sjN=0;
        trA.intronMotifs[0]=0;trA.intronMotifs[1]=0;trA.intronMotifs[2]=0;
        for (uint iex=0;iex<trA.nExons-1;iex++) {
            if (trA.canonSJ[iex]>=0)
            {//junctions - others are indels
                sjN++;
                trA.intronMotifs[trA.sjStr[iex]]++;
            };
        };

        if (trA.intronMotifs[1]>0 && trA.intronMotifs[2]==0)
            trA.sjMotifStrand=1;
        else if (trA.intronMotifs[1]==0 && trA.intronMotifs[2]>0)
            trA.sjMotifStrand=2;
        else
            trA.sjMotifStrand=0;

        if (trA.intronMotifs[1]>0 && trA.intronMotifs[2]>0 && P.outFilterIntronStrands=="RemoveInconsistentStrands")
                return;

        if (sjN>0 && trA.sjMotifStrand==0 && P.outSAMstrandField.type==1) {//strand not defined for a junction
            return;
        };

        if (P.outFilterIntronMotifs=="None") {//no filtering

        } else if (P.outFilterIntronMotifs=="RemoveNoncanonical") {
            for (uint iex=0;iex<trA.nExons-1;iex++) {
                if (trA.canonSJ[iex]==0) return;
            };
        } else if (P.outFilterIntronMotifs=="RemoveNoncanonicalUnannotated") {
            for (uint iex=0;iex<trA.nExons-1;iex++) {
                if (trA.canonSJ[iex]==0 && trA.sjAnnot[iex]==0) return;
            };
        } else {
            ostringstream errOut;
            errOut << "EXITING because of FATAL INPUT error: unrecognized value of --outFilterIntronMotifs=" <<P.outFilterIntronMotifs <<"\n";
            errOut << "SOLUTION: re-run STAR with --outFilterIntronMotifs = None -OR- RemoveNoncanonical -OR- RemoveNoncanonicalUnannotated\n";
            exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        };

        {//check mapped length for each mate
            uint nsj=0,exl=0;
            for (uint iex=0;iex<trA.nExons;iex++) {//
                exl+=trA.exons[iex].L;
                if (iex==trA.nExons-1 || trA.canonSJ[iex]==-3) {//mate is completed, make the checks
                    if (nsj>0 && (exl<P.alignSplicedMateMapLmin || exl < (uint) (P.alignSplicedMateMapLminOverLmate*RA->readLength[trA.exons[iex].iFrag])) ) {
                        return; //do not record this transcript
                    };
                    exl=0;nsj=0;
                } else if (trA.canonSJ[iex]>=0) {
                    nsj++;
                };
            };
        };

        if (P.outFilterBySJoutStage==2) {//junctions have to be present in the filtered set P.sjnovel
            for (uint iex=0;iex<trA.nExons-1;iex++) {
                if (trA.canonSJ[iex]>=0 && trA.sjAnnot[iex]==0) {
                    uint jS=trA.exons[iex].G+trA.exons[iex].L;
                    uint jE=trA.exons[iex+1].G-1;
                    if ( binarySearch2(jS,jE,P.sjNovelStart,P.sjNovelEnd,P.sjNovelN) < 0 ) return;
                };
            };
        };

        if ( trA.exons[0].iFrag!=trA.exons[trA.nExons-1].iFrag ) {//check for correct overlap between mates
            if (trA.exons[trA.nExons-1].G+trA.exons[trA.nExons-1].L <= trA.exons[0].G) return; //to avoid negative insert size
            uint iexM2=trA.nExons;
            for (uint iex=0;iex<trA.nExons-1;iex++) {//find the first exon of the second mate
                if (trA.canonSJ[iex]==-3) {//
                    iexM2=iex+1;
                    break;
                };
            };

            if ( trA.exons[iexM2-1].G + trA.exons[iexM2-1].L > trA.exons[iexM2].G ) {//mates overlap - check consistency of junctions

                if (trA.exons[0].G > \
                    trA.exons[iexM2].G+trA.exons[0].R+P.alignEndsProtrude.nBasesMax) return; //LeftMateStart > RightMateStart + allowance
                if (trA.exons[iexM2-1].G+trA.exons[iexM2-1].L > \
                   trA.exons[trA.nExons-1].G+Lread-trA.exons[trA.nExons-1].R+P.alignEndsProtrude.nBasesMax) return; //LeftMateEnd   > RightMateEnd +allowance

                //check for junctions consistency
                uint iex1=1, iex2=iexM2+1; //last exons of the junction
                for  (; iex1<iexM2; iex1++) {//find first junction that overlaps 2nd mate
                    if (trA.exons[iex1].G >= trA.exons[iex2-1].G + trA.exons[iex2-1].L) break;
                };
                while (iex1<iexM2 && iex2<trA.nExons) {//cycle through all overlapping exons
                    if (trA.canonSJ[iex1-1]<0) {//skip non-junctions
                        iex1++;
                        continue;
                    };
                    if (trA.canonSJ[iex2-1]<0) {//skip non-junctions
                        iex2++;
                        continue;
                    };

                    if ( ( trA.exons[iex1].G!=trA.exons[iex2].G ) || ( (trA.exons[iex1-1].G+trA.exons[iex1-1].L) != (trA.exons[iex2-1].G+trA.exons[iex2-1].L) ) ) {
                        return; //inconsistent junctions on overlapping mates
                    };
                    iex1++;
                    iex2++;

                };//cycle through all overlapping exons
            };//mates overlap - check consistency of junctions
        };//check for correct overlap between mates

        if (P.scoreGenomicLengthLog2scale!=0) {//add gap length score
            Score += int(ceil( log2( (double) ( trA.exons[trA.nExons-1].G+trA.exons[trA.nExons-1].L - trA.exons[0].G) ) \
                     * P.scoreGenomicLengthLog2scale - 0.5));
            Score = max(0,Score);
        };

        //calculate some final values for the transcript

        trA.roStart = (trA.roStr == 0) ? trA.rStart : Lread - trA.rStart - trA.rLength;
        trA.maxScore=Score;

        if (trA.exons[0].iFrag==trA.exons[trA.nExons-1].iFrag) {//mark single fragment transcripts
            trA.iFrag=trA.exons[0].iFrag;
            RA->maxScoreMate[trA.iFrag] = max (RA->maxScoreMate[trA.iFrag] , Score);
        } else {
            trA.iFrag=-1;
        };

        //Variation
        //Score+=trA.variationAdjust(mapGen, R);

        trA.maxScore=Score;

        // transcript has been finalized, compare the score and record
        if (       Score+P.outFilterMultimapScoreRange >= wTr[0]->maxScore \
                || ( trA.iFrag>=0 && Score+P.outFilterMultimapScoreRange >= RA->maxScoreMate[trA.iFrag] ) \
                || P.pCh.segmentMin>0) {
                //only record the transcripts within the window that are in the Score range
                //OR within the score range of each mate
                //OR all transcript if chimeric detection is activated

//             if (P.alignEndsType.in=="EndToEnd") {//check that the alignment is end-to-end
//                 uint rTotal=trA.rLength+trA.lIns;
// //                 for (uint iex=1;iex<trA.nExons;iex++) {//find the inside exons
// //                     rTotal+=trA.exons[iex].R-trA.exons[iex-1].R;
// //                 };
//                 if ( (trA.iFrag<0 && rTotal<(RA->readLength[0]+RA->readLength[1])) || (trA.iFrag>=0 && rTotal<RA->readLength[trA.iFrag])) return;
//             };

            uint iTr=0; //transcript insertion/replacement place

            trA.mappedLength=0;
            for (uint iex=0;iex<trA.nExons;iex++) {//caclulate total mapped length
                trA.mappedLength += trA.exons[iex].L;
            };

            while (iTr < *nWinTr) {//scan through all recorded transcripts for this window - check for duplicates

                //another way to calculate uOld, uNew: w/o gMap
                uint nOverlap=blocksOverlap(trA,*wTr[iTr]);
                uint uNew=trA.mappedLength-nOverlap;
                uint uOld=wTr[iTr]->mappedLength-nOverlap;

                if (uNew==0 && Score < wTr[iTr]->maxScore) {//new transript is a subset of the old ones
                    break;
                } else if (uOld==0) {//old transcript is a subset of the new one, remove old transcript
                    Transcript *pTr=wTr[iTr];
                    for  (uint ii=iTr+1;ii<*nWinTr;ii++) wTr[ii-1]=wTr[ii]; //shift transcripts
                    (*nWinTr)--;
                    wTr[*nWinTr]=pTr;
                } else if (uOld>0 && (uNew>0 || Score >= wTr[iTr]->maxScore) ) {//check next transcript
                    iTr++;
                };

            };

            if (iTr==*nWinTr) {//insert the new transcript
                for (iTr=0;iTr<*nWinTr;iTr++) {//find inseriton location
                    if (Score>wTr[iTr]->maxScore || (Score==wTr[iTr]->maxScore && trA.gLength<wTr[iTr]->gLength) ) break;
                };

                Transcript *pTr=wTr[*nWinTr];
                for (int ii=*nWinTr; ii> int(iTr); ii--) {//shift all the transcript pointers down from iTr
                    wTr[ii]=wTr[ii-1];
                };
                wTr[iTr]=pTr; //the new transcript pointer is now at *nWinTr+1, move it into the iTr
                *(wTr[iTr])=trA;
                if (*nWinTr<P.alignTranscriptsPerWindowNmax) {
                    (*nWinTr)++; //increment number of transcripts per window;
                } else {
                        //"WARNING: too many recorded transcripts per window: iRead="<<RA->iRead<< "\n";
                };
            };
        };


        return;
    };

    ///////////////////////////////////////////////////////////////////////////////////
    int dScore=0;
    Transcript trAi=trA; //trA copy with this align included, to be used in the 1st recursive call of StitchAlign
    if (trA.nExons>0) {//stitch, a transcript has already been originated

        dScore=stitchAlignToTranscript(tR2, tG2, WA[iA].rStart, WA[iA].gStart, WA[iA].Length, WA[iA].iFrag,  WA[iA].sjA, P, R, mapGen, &trAi, RA->outFilterMismatchNmaxTotal);
        //TODO check if the new stitching creates too many MM, quit this transcript if so

    } else { //this is the first align in the transcript
            trAi.exons[0].R=trAi.rStart=WA[iA].rStart; //transcript start/end
            trAi.exons[0].G=trAi.gStart=WA[iA].gStart;
            trAi.exons[0].L=WA[iA].Length;
            trAi.exons[0].iFrag=WA[iA].iFrag;
            trAi.exons[0].sjA=WA[iA].sjA;

            trAi.nExons=1; //recorded first exon

            for (uint ii=0;ii<WA[iA].Length;ii++) dScore+=scoreMatch; //sum all the scores

            trAi.nMatch=WA[iA].Length; //# of matches

            for (uint ii=0; ii<nA; ii++) WAincl[ii]=false;


    };

    if (dScore>-1000000) {//include this align
        WAincl[iA]=true;

        if ( WA[iA].Nrep==1 ) trAi.nUnique++; //unique piece
        if ( WA[iA].Anchor>0 ) trAi.nAnchor++; //anchor piece

        stitchWindowAligns(iA+1, Score+dScore, WAincl, WA[iA].rStart+WA[iA].Length-1, WA[iA].gStart+WA[iA].Length-1, trAi, Lread, WA, R, mapGen, P, wTr, nWinTr, RA);
    } else {

    };

    //also run a transcript w/o including this align
    if (WA[iA].Anchor!=2 || trA.nAnchor>0) {//only allow exclusion if this is not the last anchor, or other anchors have been used
        WAincl[iA]=false;
        stitchWindowAligns(iA+1, Score, WAincl, tR2, tG2, trA, Lread, WA, R, mapGen, P, wTr, nWinTr, RA);
    };
    return;
};


