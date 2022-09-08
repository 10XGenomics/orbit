#include "IncludeDefine.h"
#include "Parameters.h"
#include "ReadAlign.h"
#include "ErrorWarning.h"

void ReadAlign::assignAlignToWindow(uint a1, uint aLength, uint aStr, uint aNrep, uint aFrag, uint aRstart,
                                    bool aAnchor, uint sjA) {

    uint iW=winBin[aStr][a1>>P.winBinNbits];

    if (iW==uintWinBinMax || (!aAnchor && aLength < WALrec[iW]) ) return; //alignment does not belong to any window, or it's shorter than rec-length

    //check if this alignment overlaps with any other alignment in the window, record the longest of the two
    {//do not check for overlap if this is an sj-align
        uint iA;
        for (iA=0; iA<WA[iW].size(); iA++) {
            if (aFrag==WA[iW][iA].iFrag && WA[iW][iA].sjA==sjA \
                && a1+WA[iW][iA].rStart==WA[iW][iA].gStart+aRstart \
                && ( (aRstart>=WA[iW][iA].rStart && aRstart<WA[iW][iA].rStart+WA[iW][iA].Length) \
                  || (aRstart+aLength>=WA[iW][iA].rStart && aRstart+aLength<WA[iW][iA].rStart+WA[iW][iA].Length) ) ) {//this piece overlaps with iA
                break;
            };
        };
        if (iA<WA[iW].size()) {//found overlap
            if (aLength>WA[iW][iA].Length) {//replace

                uint iA0;//iA0 is where the align has to be inserted
                for (iA0=0;iA0<WA[iW].size();iA0++)
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
                        WA[iW][iA1]=WA[iW][iA1-1];
                    };
                } else if (iA0>iA) {//shift aligns up
                    for (uint iA1=iA;iA1<iA0;iA1++) {//shift aligns to free up insertion point
                        WA[iW][iA1]=WA[iW][iA1+1];
                    };
                };

                WA[iW][iA0] = uiWA(aLength, aRstart, a1,
                                   aNrep, int(aAnchor), aFrag, sjA);

            };
            return; //do not record new align
        };
    };


    if (WA[iW].size()==P.seedPerWindowNmax) {//too many aligns per window,  re-calcualte min-length, remove the shortest one,

        WALrec[iW]=Lread+1;
        for (uint iA=0; iA<WA[iW].size(); iA++) {//find the new min-length
            if (WA[iW][iA].Anchor!=1) WALrec[iW]=min(WALrec[iW],WA[iW][iA].Length); //protect the anchors - they are not counted for min-length
        };


        if (WALrec[iW]==Lread+1) {//this could happen if there are too many anchors
            mapMarker=MARKER_TOO_MANY_ANCHORS_PER_WINDOW;
            WC.clear();
            return;
        };


        if (!aAnchor && aLength < WALrec[iW]) return; //alignment is shorter than min-length, do not record - unless it's an anchor

        uint iA1=0;
        for (uint iA=0; iA<WA[iW].size(); iA++) {//remove the shortest aligns
            if ( WA[iW][iA].Anchor==1 || WA[iW][iA].Length > WALrec[iW] ) {//re-record the anchors and long aligns
                WA[iW][iA1]=WA[iW][iA]; //re-record the iA-th alignment into iA1-th place
                iA1++;
            };
        };
        WA[iW].resize(iA1);

        if (!aAnchor && aLength <= WALrec[iW]) {//current align was removed, zero out its nWAP
            nWAP[iW]=0;
        };

    };
    auto& cWA = WA[iW];
    if ( aAnchor || aLength > WALrec[iW] ) {
        if (cWA.size() >= P.seedPerWindowNmax) {
            exitWithError("BUG: iA>=P.seedPerWindowNmax in stitchPieces, exiting",std::cerr, P.inOut->logMain, EXIT_CODE_BUG, P);
        };

        uint iA;
        for (iA = 0; iA < cWA.size(); iA++) {//find the insertion point in case aligns are not sorted by aRstart
            //TODO binary search
            if (aRstart < cWA[iA].rStart) break;
        };
        cWA.insert(cWA.begin() + iA1, 1,
                   uiWA(aLength, aRstart, a1,
                        aNrep, int(aAnchor), aFrag, sjA));
        nWAP[iW]++;
        if (aAnchor && WlastAnchor[iW]<iA) WlastAnchor[iW]=iA; //record the index of the last anchor
    };
};
