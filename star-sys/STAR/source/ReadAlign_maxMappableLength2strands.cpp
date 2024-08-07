#include <array>
#include <vector>

#include "ReadAlign.h"
#include "SuffixArrayFuns.h"
#include "ErrorWarning.h"

using std::array;
using std::vector;

uint ReadAlign::maxMappableLength2strands(uint pieceStartIn, uint pieceLengthIn, uint iDir, uint iSA1, uint iSA2, uint& maxLbest, uint iFrag) {
    //returns number of mappings, maxMappedLength=mapped length
    uint Nrep=0, indStartEnd[2], maxL;

    vector<uint> NrepAll(P.pGe.gSAsparseD);
    vector<array<uint, 2>> indStartEndAll(P.pGe.gSAsparseD);
    vector<uint> maxLall(P.pGe.gSAsparseD);
    maxLbest=0;

    bool dirR = iDir==0;

    // defaults:  (from genomeParameters.txt)
    // gSAsparseD = 1
    // gSAindexNbases = 14

    for (uint iDist=0; iDist<min(pieceLengthIn,P.pGe.gSAsparseD); iDist++) {//cycle through different distances
        uint pieceStart;
        uint pieceLength=pieceLengthIn-iDist;

        //calculate full index
        uint Lmax=min(P.pGe.gSAindexNbases,pieceLength);
        uint ind1=0;
        if (dirR) {//forward search
            pieceStart=pieceStartIn+iDist;
            for (uint ii=0;ii<Lmax;ii++) {//calculate index TODO: make the index calculation once for the whole read and store it
                ind1 <<=2LLU;
                ind1 += ((uint) Read1[0][pieceStart+ii]);
            };
        } else {//reverse search
            pieceStart=pieceStartIn-iDist;
            for (uint ii=0;ii<Lmax;ii++) {//calculate index TODO: make the index calculation once for the whole read and store it
                ind1 <<=2LLU;
                ind1 += ( 3-((uint) Read1[0][pieceStart-ii]) );
            };
        };

        //find SA boundaries
        uint Lind=Lmax;
        while (Lind>0) {//check the presence of the prefix for Lind
            iSA1=mapGen.SAi[mapGen.genomeSAindexStart[Lind-1]+ind1]; // starting point for suffix array search.
            if ((iSA1 & mapGen.SAiMarkAbsentMaskC) == 0) {//prefix exists
                break;
            } else {//this prefix does not exist, reduce Lind
                --Lind;
                ind1 = ind1 >> 2;
            };
        };

        // define upper bound for suffix array range search.
        bool iSA2good = true;
        if (mapGen.genomeSAindexStart[Lind-1]+ind1+1 < mapGen.genomeSAindexStart[Lind]) {//we are not at the end of the SA
            iSA2 = mapGen.SAi[mapGen.genomeSAindexStart[Lind-1]+ind1+1];
            if ( (iSA2 & mapGen.SAiMarkAbsentMaskC) == 0) {
                iSA2 = (iSA2 & mapGen.SAiMarkNmask) - 1;
            } else {
                iSA2 = mapGen.nSA-1; //safe, but can probably do better
                iSA2good = false;
            };
        } else {
            iSA2=mapGen.nSA-1;
            iSA2good = false;
        };

    //#define SA_SEARCH_FULL

    #ifdef SA_SEARCH_FULL
        //full search of the array even if the index search gave maxL
        maxL=0;
        Nrep = maxMappableLength(mapGen, Read1, pieceStart, pieceLength, iSA1 & mapGen.SAiMarkNmask, iSA2, dirR, maxL, indStartEnd);
    #else
        bool iSA1noN = (iSA1 & mapGen.SAiMarkNmaskC)==0;
        if (Lind < P.pGe.gSAindexNbases && iSA1noN && iSA2good) {//no need for SA search
            // very short seq, already found hits in suffix array w/o having to search the genome for extensions.
            indStartEnd[0]=iSA1;
            indStartEnd[1]=iSA2;
            Nrep=indStartEnd[1]-indStartEnd[0]+1;
            maxL=Lind;
        } else if (iSA1==iSA2 && iSA1noN && iSA2good) {//unique align already, just find maxL
            if ((iSA1 & mapGen.SAiMarkNmaskC)!=0) {
                ostringstream errOut;
                errOut  << "BUG: in ReadAlign::maxMappableLength2strands";
                exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_BUG, P);
            };
            indStartEnd[0]=indStartEnd[1]=iSA1;
            Nrep=1;
            bool comparRes;
            maxL=compareSeqToGenome(mapGen, Read1.data(), pieceStart, pieceLength, Lind, iSA1, dirR, comparRes);
        } else {//need SA search, pieceLength>maxL
            if (iSA2good && iSA1noN) {
                maxL = Lind; //Lind bases were already matched
            } else {
                maxL=0;
            };        
            Nrep = maxMappableLength(mapGen, Read1.data(), pieceStart, pieceLength, iSA1 & mapGen.SAiMarkNmask, iSA2, dirR, maxL, indStartEnd);
        };
    #endif

        if (maxL+iDist > maxLbest) {//this idist is better
            maxLbest=maxL+iDist;
        };
        NrepAll[iDist]=Nrep;
        indStartEndAll[iDist][0]=indStartEnd[0];
        indStartEndAll[iDist][1]=indStartEnd[1];
        maxLall[iDist]=maxL;
    };

    for (uint iDist=0; iDist<min(pieceLengthIn,P.pGe.gSAsparseD); iDist++) {//cycle through different distances, store the ones with largest maxL
        if ( (maxLall[iDist]+iDist) == maxLbest) {
            storeAligns(iDir, (dirR ? pieceStartIn+iDist : pieceStartIn-iDist), NrepAll[iDist], maxLall[iDist], indStartEndAll[iDist].data(), iFrag);
        };
    };
    return Nrep;
};
