#include "Transcriptome.h"
#include "ReadAlign.h"
#include "serviceFuns.cpp"


int alignToTranscript(Transcript &aG, uint trS1, uint8 trStr1, uint32 *exSE1, uint32 *exLenCum1, uint16 exN1, Transcript &aT) {

    //find exon that overlaps beginning of the read
    uint32 g1=aG.exons[0].G-trS1;//start of the transcript
    uint32 ex1=binarySearch1<uint32>(g1, exSE1, 2*exN1);
    if (ex1>=2*exN1) return 0; //align start is to the right of all exons

    if (ex1%2==1) {//beginning of the read >=end of an exon
        if (exSE1[ex1]==g1) {//first base of the read is exactly the last base of the exon
            --ex1;
        } else {
            return 0;//beginning of the read is past the end of an exon, align does not belong to this transcript
        };
    };
    ex1=ex1/2; //this is the first exon of the alignment

    aT.nExons=0;
    aT.primaryFlag=false;

    aG.canonSJ[aG.nExons-1]=-999; //marks the last exons
    for (uint32 iab=0; iab<aG.nExons; iab++) {//scan through all blocks of the align
        if (aG.exons[iab].G+aG.exons[iab].L>exSE1[2*ex1+1]+trS1+1) {//block extends past exon end
            return 0;
        };

        if (iab==0 || aG.canonSJ[iab-1]<0) {
            aT.exons[aT.nExons].R=aG.exons[iab].R;
            aT.exons[aT.nExons].G=aG.exons[iab].G-trS1-exSE1[2*ex1]+exLenCum1[ex1];
            aT.exons[aT.nExons].L=aG.exons[iab].L;
            aT.exons[aT.nExons].iFrag=aG.exons[iab].iFrag;
            if (aT.nExons>0) aT.canonSJ[aT.nExons-1]=aG.canonSJ[iab-1];
            ++aT.nExons;
        } else {
            aT.exons[aT.nExons-1].L+=aG.exons[iab].L;
        };
        switch (aG.canonSJ[iab]) {
            case -999: //last exon
                if (trStr1==2) {//convert align coordinates if on the -strand
                    uint32 trlength=exLenCum1[exN1-1]+exSE1[2*exN1-1]-exSE1[2*exN1-2]+1; //transcript length
                    for (uint32 iex=0; iex<aT.nExons; iex++) {
                        aT.exons[iex].R=aG.Lread-(aT.exons[iex].R+aT.exons[iex].L);
                        aT.exons[iex].G=trlength-(aT.exons[iex].G+aT.exons[iex].L);  // slow loop
                    };
                    for (uint32 iex=0; iex<aT.nExons/2; iex++) {
                        swap(aT.exons[iex].R,aT.exons[aT.nExons-1-iex].R); // also slow
                        swap(aT.exons[iex].G,aT.exons[aT.nExons-1-iex].G);
                        swap(aT.exons[iex].L,aT.exons[aT.nExons-1-iex].L);
                        swap(aT.exons[iex].iFrag,aT.exons[aT.nExons-1-iex].iFrag);
                    };
                    for (uint32 iex=0; iex<(aT.nExons-1)/2; iex++) {
                        swap(aT.canonSJ[iex],aT.canonSJ[aT.nExons-2-iex]);
                    };
                };
                for (uint32 iex=0; iex<aT.nExons; iex++) {//no junctions in the transcritomic coordinates
                    aT.sjAnnot[iex]=0;
                    aT.shiftSJ[iex][0]=0;
                    aT.shiftSJ[iex][1]=0;
                    aT.sjStr[iex]=0;
                };

                return 1; //reached the end of blocks, align is consistent with this transcript
                break;
            case -3: //mate connection
                ex1=binarySearch1<uint32>(aG.exons[iab+1].G-trS1, exSE1, 2*exN1);
                if (ex1%2==1) {//beginning of the mext mate in the middle of the exon?
                    return 0; //align does not belong to this transcript
                } else {
                    ex1=ex1/2; //this is the first exon of the second mate
                };
                break;
            case -2: //insertion
                break;
            case -1: //deletion
                break;
            default://junctions
                if ( aG.exons[iab].G+aG.exons[iab].L==exSE1[2*ex1+1]+trS1+1 && aG.exons[iab+1].G==exSE1[2*(ex1+1)]+trS1 ) {
                    //junction matches transcript junction
                    ++ex1;
                } else {
                    return 0;
                };
        };
    };
    return 0; //this should not happen
};

uint32 Transcriptome::quantAlign (Transcript &aG, Transcript *aTall, vector<uint32> &/*readTranscripts*/, set<uint32> &/*readTrGenes*/) {
    uint32 nAtr=0; //number of alignments to the transcriptome

    //binary search through transcript starts
    uint32 tr1=binarySearch1a<uint>(aG.exons[0].G, trS, nTr);
    if (tr1==(uint32) -1) return 0; //alignment outside of range of all transcripts

    uint aGend=aG.exons[aG.nExons-1].G;

    ++tr1;
    do {//cycle back through all the transcripts
        --tr1;
        if (aGend<=trE[tr1]) {//this transcript contains the read
                int aStatus=alignToTranscript(aG, trS[tr1], trStr[tr1], exSE+2*trExI[tr1], exLenCum+trExI[tr1], trExN[tr1], aTall[nAtr]);
                if (aStatus==1) {//align conforms with the transcript
                    aTall[nAtr].Chr = tr1;
                    aTall[nAtr].Str = trStr[tr1]==1 ? aG.Str : 1-aG.Str; //TODO strandedness
                    ++nAtr;
                };
        };
    } while (trEmax[tr1]>=aGend && tr1>0);

    return nAtr;
};
