#include "IncludeDefine.h"
#include "ReadAlign.h"
#include "stitchWindowAligns.h"
#include "sjAlignSplit.cpp"
#include "alignSmithWaterman.h"

void ReadAlign::stitchPieces(char **R, uint Lread) {

    //zero-out winBin
    memset(winBin[0],255,sizeof(winBin[0][0])*P.winBinN);
    memset(winBin[1],255,sizeof(winBin[0][0])*P.winBinN);

//     for (uint iWin=0;iWin<nWall;iWin++) {//zero out winBin
//         if (WC[iWin][WC_gStart]<=WC[iWin][WC_gEnd]) {//otherwise the window is dead
//             memset(&(winBin[WC[iWin][WC_Str]][WC[iWin][WC_gStart]]),255,sizeof(winBin[0][0])*(WC[iWin][WC_gEnd]-WC[iWin][WC_gStart]+1));
//         };
// //         for (uint ii=C[iWin][WC_gStart]; ii<WC[iWin][WC_gEnd]; ii++) {
// //             winBin[WC[WC_Str]
// //         };
//     };

//     //debug
//     for (uint ii=0;ii<P.winBinN;ii++){
//         if (winBin[0][ii]!=uintWinBinMax || winBin[1][ii]!=uintWinBinMax) {
//             cerr<< "BUG in stitchPieces: ii="<<ii<<"   "<< winBin[0][ii] <<"   "<<winBin[1][ii] <<"   iRead="<<iRead<<"   nW="<<nW<<endl;
//             for (uint iWin=0;iWin<nW;iWin++) {
//                 cerr <<WC[iWin][WC_gStart]<<"   " <<WC[iWin][WC_gEnd] <<"   "<<WC[iWin][WC_Str] <<endl;
//             };
//             exit(1);
//         };
//     };

    WC.clear(); //number of windows
//    for (uint iP=0; iP<nP; iP++) {
//        cout << "NREP: " << PC[iP][PC_Nrep] << endl;
//        cout << "Dir: " << PC[iP][PC_Dir] << endl;
//        cout << "Length: " << PC[iP][PC_Length] << endl;
//        cout << endl;
//    }
    // Alignment windows are binned regions of the genome that anchor pieces lie in.
    for (uint iP=0; iP<nP; iP++) {//scan through all anchor pieces, create alignment windows
        // np is number of pieces (stored seed alignments)


//          if (PC[iP][PC_Nrep]<=P.winAnchorMultimapNmax || PC[iP][PC_Length]>=readLength[PC[iP][PC_iFrag]] ) {//proceed if piece is an anchor, i.e. maps few times or is long enough
       if (PC[iP][PC_Nrep]<=P.winAnchorMultimapNmax ) {//proceed if piece is an anchor, i.e. maps few times

            uint aDir   = PC[iP][PC_Dir];
            uint aLength= PC[iP][PC_Length];

            for (uint iSA=PC[iP][PC_SAstart]; iSA<=PC[iP][PC_SAend]; iSA++) {//scan through all alignments of this piece
                // going through ordered positions in the suffix array from PC_SAstart to PC_SAend
                uint a1 = mapGen.SA[iSA];
                //printf("a1 %llu\n", a1);
                uint aStr = a1 >> mapGen.GstrandBit;
                a1 &= mapGen.GstrandMask; //remove strand bit

                //convert to positive strand
                if (aDir==1 && aStr==0) {
                    aStr=1;
                } else if (aDir==0 && aStr==1) {
                    a1 = mapGen.nGenome - (aLength+a1);
                } else if (aDir==1 && aStr==1) {
                    aStr=0;
                    a1 = mapGen.nGenome - (aLength+a1);
                };
                //final strand
                if (revertStrand) { //modified strand according to user input CHECK!!!!
                    aStr=1-aStr;
                };


                if (a1>=mapGen.sjGstart) {//this is sj align
                    uint a1D, aLengthD, a1A, aLengthA, sj1;
                    if (sjAlignSplit(a1, aLength, mapGen, a1D, aLengthD, a1A, aLengthA, sj1)) {//align crosses the junction

                        int addStatus=createExtendWindowsWithAlign(a1D, aStr);//add donor piece
                        if (addStatus==EXIT_createExtendWindowsWithAlign_TOO_MANY_WINDOWS) {//too many windows
                            break;
                        };
                        addStatus=createExtendWindowsWithAlign(a1A, aStr);//add acceptor piece
                        if (addStatus==EXIT_createExtendWindowsWithAlign_TOO_MANY_WINDOWS) {//too many windows
                            break;
                        };
                    };
                } else {//this is a normal genomic read
                    int addStatus=createExtendWindowsWithAlign(a1, aStr);
                    if (addStatus==EXIT_createExtendWindowsWithAlign_TOO_MANY_WINDOWS) {//too many windows
                        break;
                    };
                };
            }; //for (uint iSA=PC[iP][PC_SAstart]; iSA<=PC[iP][PC_SAend]; iSA++) //scan through all alignments of this piece
        };//if (PC[iP][PC_Nrep]<=P.winAnchorMultimapNmax) //proceed if anchor
    };//for (uint iP=0; iP<nP; iP++) //scan through all anchor pieces, create alignment windows


    // Pass through all windows and add flanking regions
    for (uint iWin=0;iWin<WC.size();iWin++) {//extend windows with flanks
        if (WC[iWin].gStart<=WC[iWin].gEnd) {//otherwise the window is dead, happens when it is merged to a neighboring window

            uint wb=WC[iWin].gStart;
            for (uint ii=0; ii<P.winFlankNbins && wb>0 && mapGen.chrBin[(wb-1) >> P.winBinChrNbits]==WC[iWin].Chr;ii++) {
                wb--;
                winBin[ WC[iWin].Str ][ wb ]=(uintWinBin) iWin;
            };
            WC[iWin].gStart = wb;

            wb=WC[iWin].gEnd;
            for (uint ii=0; ii<P.winFlankNbins && wb+1<P.winBinN && mapGen.chrBin[(wb+1) >> P.winBinChrNbits]==WC[iWin].Chr;ii++) {
                wb++;
                winBin[ WC[iWin].Str ][ wb ]=(uintWinBin) iWin;
            };
            WC[iWin].gEnd = wb;


        };
        WA[iWin].clear(); //initialize nWA
        WALrec[iWin]=0; //initialize rec-length
        WlastAnchor[iWin]=-1;
    };

    #ifdef OFF_BEFORE_SEEDdistribution
        #warning OFF_BEFORE_SEEDdistribution
        nW=0;
        nTr=0;
        return;
    #endif

    //scan through all pieces/aligns, add them to alignment windows, create alignment coordinates
    for (uint iP=0; iP<nP; iP++) {
        uint aNrep=PC[iP][PC_Nrep];
        uint aFrag=PC[iP][PC_iFrag];
        uint aLength=PC[iP][PC_Length];
        uint aDir=PC[iP][PC_Dir];

        bool aAnchor=(aNrep<=P.winAnchorMultimapNmax); //this align is an anchor or not

        for (uint ii=0;ii<WC.size();ii++) {//initialize nWAP (number per window piece?)
            nWAP[ii]=0;
        };


        for (uint iSA=PC[iP][PC_SAstart]; iSA<=PC[iP][PC_SAend]; iSA++) {//scan through all alignments

            uint a1 = mapGen.SA[iSA];
            uint aStr = a1 >> mapGen.GstrandBit;
            a1 &= mapGen.GstrandMask; //remove strand bit
            uint aRstart=PC[iP][PC_rStart];

            //convert to positive strand
            if (aDir==1 && aStr==0) {
                aStr=1;
                aRstart = Lread - (aLength+aRstart);
            } else if (aDir==0 && aStr==1) {
                aRstart = Lread - (aLength+aRstart);
                a1 = mapGen.nGenome - (aLength+a1);
            } else if (aDir==1 && aStr==1) {
                aStr=0;
                a1 = mapGen.nGenome - (aLength+a1);
            };

            //final strand
            if (revertStrand) { //modified strand according to user input CHECK!!!!
                aStr=1-aStr;
            };


            if (a1>=mapGen.sjGstart) {//this is sj read
                uint a1D, aLengthD, a1A, aLengthA, isj1;
                if (sjAlignSplit(a1, aLength, mapGen, a1D, aLengthD, a1A, aLengthA, isj1)) {//align crosses the junction

                        assignAlignToWindow(a1D, aLengthD, aStr, aNrep, aFrag, aRstart, aAnchor, isj1);
                        assignAlignToWindow(a1A, aLengthA, aStr, aNrep, aFrag, aRstart+aLengthD, aAnchor, isj1);

                  } else {//align does not cross the junction
                        continue; //do not check this align, continue to the next one
                  };

              } else {//this is a normal genomic read
                    assignAlignToWindow(a1, aLength, aStr, aNrep, aFrag, aRstart, aAnchor, -1);
              };
        };


//         for (uint ii=0;ii<nW;ii++) {//check of some pieces created too many aligns in some windows, and remove those from WA (ie shift nWA indices
//             if (nWAP[ii]>P.seedNoneLociPerWindow) nWA[ii] -= nWAP[ii];
//         };
    };

    //generate transcript for each window, choose the best
    trBest =&trInit; //initialize next/best
    uint iW1=0;//index of non-empty windows
    uint trNtotal=0; //total number of recorded transcripts


    for (uint iW=0; iW<WC.size(); iW++) {//transcripts for all windows

        if (WA[iW].size()==0) continue; //the window does not contain any aligns because it was merged with other windows

        if (WlastAnchor[iW]<WA[iW].size()) {
            WA[ iW ][ WlastAnchor[iW] ].Anchor=2; //mark the last anchor
        };
        // TODO: The `.include` value likely is only relevant for a WA[iW] row, consider
        //  making that more specific and not allocating a bool for all elements.
        for (uint ii=0;ii<WA[iW].size();ii++) WA[iW][ii].include=false; //initialize mask

        trA=trInit; //that one is initialized
        trA.Chr = WC[iW].Chr;
        trA.Str = WC[iW].Str;
        trA.roStr = revertStrand ? 1-trA.Str : trA.Str; //original strand of the read
        trA.maxScore=0;
        trAll[iW1]=trArrayPointer+trNtotal;
        if (trNtotal+P.alignTranscriptsPerWindowNmax >= P.alignTranscriptsPerReadNmax) {
            P.inOut->logMain << "WARNING: not enough space allocated for transcript. Did not process all windows for read "<< readName+1 <<endl;
            P.inOut->logMain <<"   SOLUTION: increase alignTranscriptsPerReadNmax and re-run\n" << flush;
            break;
        };
        //printf("trA %llu\n", trA.Chr);
        *(trAll[iW1][0])=trA;
        nWinTr[iW1]=0; //initialize number of transcripts per window

        stitchWindowAligns(0, WA[iW].size(), 0, 0, 0, trA, Lread, WA[iW], R[trA.roStr==0 ? 0:2], mapGen, P, trAll[iW1], nWinTr+iW1, this);

        if (nWinTr[iW1]==0) {
            continue;
        };

        if (trAll[iW1][0]->maxScore > trBest->maxScore || (trAll[iW1][0]->maxScore == trBest->maxScore && trAll[iW1][0]->gLength < trBest->gLength ) ) {
            trBest=trAll[iW1][0];
        };
        //printf("why not add here %llu\n", nWinTr[iW1]);
        trNtotal += nWinTr[iW1];
        iW1++;
    };

    WC.resize(iW1);//only count windows that had alignments

//     {//debug
//         std::time(&timeFinish);
//         double timeDiff=difftime(timeFinish,timeStart);
//         cout << "     "<< timeDiff << "     "<<trBest->maxScore*100/Lread<<"   "<<iRead<<endl;;
//     };

    if (trBest->maxScore==0) {//no window was aligned (could happen if for all windows too many reads are multiples)
        mapMarker = MARKER_NO_GOOD_WINDOW;
        WC.clear();
        return;
    };

};//end of function
