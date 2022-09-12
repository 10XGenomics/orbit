#include "IncludeDefine.h"
#include "Parameters.h"
#include "ReadAlign.h"
#include "ErrorWarning.h"

void ReadAlign::assignAlignToWindow(uint a1, uint aLength, uint aStr, uint aNrep, uint aFrag, uint aRstart,
                                    bool aAnchor, uint sjA) {

    uint iW=winBin[aStr][a1>>P.winBinNbits];

    if (iW==uintWinBinMax || (!aAnchor && aLength < WALrec[iW]) ) return; //alignment does not belong to any window, or it's shorter than rec-length
    auto& cWA = WA[iW];
    //check if this alignment overlaps with any other alignment in the window, record the longest of the two
    {//do not check for overlap if this is an sj-align
        uint iA;
        for (iA=0; iA<cWA.size(); iA++) {
            if (aFrag==cWA[iA].iFrag && cWA[iA].sjA==sjA \
            // The new alignment genome start + old alignment read start = old alignment genome start + new read start
            // indicates one is a suffix or prefix of the other in the same position
                && a1+cWA[iA].rStart==cWA[iA].gStart+aRstart \
                //this piece overlaps with iA
                && (aRstart<=cWA[iA].rStart+cWA[iA].Length && cWA[iA].rStart <= aRstart+aLength)) {
                break;
            };
        };
        if (iA<cWA.size()) {//found overlap
            if (aLength>cWA[iA].Length) {//replace

                uint iA0;//iA0 is where the align has to be inserted
                for (iA0=0;iA0<cWA.size();iA0++)
                {//find the insertion point TODO binary search
                    if (iA0!=iA && aRstart<WA[iW][iA0].rStart)
                    {//do not compare with the piece to be removed
                        break;
                    };
                };

                if (iA0>iA)
                {//true insertion place since iA will be removed
                    --iA0;
                };

                if (iA0<iA) {//shift aligns down
                    for (uint iA1=iA;iA1>iA0;iA1--) {//shift aligns to free up insertion point
                        cWA[iA1]=cWA[iA1-1];
                    };
                } else if (iA0>iA) {//shift aligns up
                    for (uint iA1=iA;iA1<iA0;iA1++) {//shift aligns to free up insertion point
                        cWA[iA1]=cWA[iA1+1];
                    };
                };

                cWA[iA0] = uiWA(aLength, aRstart, a1,
                                   aNrep, int(aAnchor), aFrag, sjA);
            };
            return; //do not record new align
        };
    };


    if (cWA.size()==P.seedPerWindowNmax) {//too many aligns per window,  re-calcualte min-length, remove the shortest one,

        WALrec[iW]=Lread+1;
        for (uint iA=0; iA<cWA.size(); iA++) {//find the new min-length
            if (cWA[iA].Anchor!=1) WALrec[iW]=min(WALrec[iW],cWA[iA].Length); //protect the anchors - they are not counted for min-length
        };


        if (WALrec[iW]==Lread+1) {//this could happen if there are too many anchors
            mapMarker=MARKER_TOO_MANY_ANCHORS_PER_WINDOW;
            WC.clear();
            return;
        };


        if (!aAnchor && aLength < WALrec[iW]) return; //alignment is shorter than min-length, do not record - unless it's an anchor

        uint iA1=0;
        for (uint iA=0; iA<cWA.size(); iA++) {//remove the shortest aligns
            if ( cWA[iA].Anchor==1 || cWA[iA].Length > WALrec[iW] ) {//re-record the anchors and long aligns
                cWA[iA1]=cWA[iA]; //re-record the iA-th alignment into iA1-th place
                iA1++;
            };
        };
        cWA.resize(iA1);

        if (!aAnchor && aLength <= WALrec[iW]) {//current align was removed, zero out its nWAP
            nWAP[iW]=0;
        };

    };
    if ( aAnchor || aLength > WALrec[iW] ) {
        if (cWA.size() >= P.seedPerWindowNmax) {
            exitWithError("BUG: iA>=P.seedPerWindowNmax in stitchPieces, exiting",std::cerr, P.inOut->logMain, EXIT_CODE_BUG, P);
        };

        uint iA;
        for (iA = 0; iA < cWA.size(); iA++) {//find the insertion point in case aligns are not sorted by aRstart
            //TODO binary search
            if (aRstart < cWA[iA].rStart) break;
        };
        cWA.insert(cWA.begin() + iA, 1,
                   uiWA(aLength, aRstart, a1,
                        aNrep, int(aAnchor), aFrag, sjA));
        nWAP[iW]++;
        if (aAnchor && WlastAnchor[iW]<iA) WlastAnchor[iW]=iA; //record the index of the last anchor
    };
};
