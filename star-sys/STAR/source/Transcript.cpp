#include "Transcript.h"

void Transcript::reset() {
    extendL=0;
    primaryFlag=false;

    rStart=0; roStart=0; rLength=0; gStart=0; gLength=0; //read and genomic coordinates

    maxScore=0;
    nMatch=0;
    nMM=0;

    nGap=0; lGap=0; lDel=0; lIns=0; nDel=0; nIns=0;

    nUnique=nAnchor=0;
}

void Transcript::add(const Transcript *trIn) noexcept {
    maxScore+=trIn->maxScore;
    nMatch+=trIn->nMatch;
    nMM+=trIn->nMM;
    nGap+=trIn->nGap; lGap+=trIn->lGap;
    lDel+=trIn->lDel; nDel+=trIn->nDel;
    lIns+=trIn->lIns; nIns+=trIn->nIns;
    nUnique+=trIn->nUnique;
}

bool Transcript::extractSpliceJunctions(vector<array<uint64,2>> &sjOut) const
{
    bool annotYes=true;
    for (uint64 iex=0; iex<nExons-1; iex++) {//record all junctions
        if (canonSJ[iex]>=0) {//only record junctions, not indels or mate gap
            array<uint64,2> sj;
            sj[0]=exons[iex][EX_G]+exons[iex][EX_L];//start
            sj[1]=exons[iex+1][EX_G] - sj[0]; //gap
            sjOut.push_back(sj);
            if (sjAnnot[iex]==0)
                annotYes=false;//if one of the SJs is unannoated, annotYes=false
        };
    };
    return annotYes;
}

