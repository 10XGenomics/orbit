#include <vector>
#include "IncludeDefine.h"
#include "Parameters.h"
#include "Transcript.h"
#include "ReadAlign.h"

ReadAlign::ReadAlign (const Parameters& Pin, const Genome &genomeIn, Transcriptome *TrIn, int iChunk)
                    : mapGen(genomeIn), readFastq{nullptr, nullptr}, P(Pin), chunkTr(TrIn), trInit()
{
    iRead = 0;
    readFilesIndex = 0;
    readNmates=P.readNmates;
    winBin = new uintWinBin* [2];
    winBin[0] = new uintWinBin [P.winBinN];
    winBin[1] = new uintWinBin [P.winBinN];
    memset(winBin[0],255,sizeof(winBin[0][0])*P.winBinN);
    memset(winBin[1],255,sizeof(winBin[0][0])*P.winBinN);
    //RNGs
    rngMultOrder.seed(P.runRNGseed*(iChunk+1));
    rngUniformReal0to1=std::uniform_real_distribution<double> (0.0, 1.0);
    //transcriptome
    if ( P.quant.trSAM.yes ) {
        alignTrAll=new Transcript [P.alignTranscriptsPerReadNmax];
    };
    //split
    splitR=new uint*[3];
    splitR[0]=new uint[P.maxNsplit]; splitR[1]=new uint[P.maxNsplit]; splitR[2]=new uint[P.maxNsplit];
    //alignments
    PC=new uiPC[P.seedPerReadNmax];
    WC= std::vector<Window>(P.alignWindowsPerReadNmax);
    nWAP=new uint[P.alignWindowsPerReadNmax];
    WALrec=new uint[P.alignWindowsPerReadNmax];
    WlastAnchor=new uint[P.alignWindowsPerReadNmax];

    // Create a P.alignWindowsPerReadNmax x P.seedPerWindowNmax by reserving but not filling all the memory
    WA = vector<vector<uiWA>>(P.alignWindowsPerReadNmax);
    std::for_each(WA.begin(), WA.end(), [this](std::vector<uiWA> &row) { row.reserve(P.seedPerWindowNmax); });
    trAll = new Transcript**[P.alignWindowsPerReadNmax+1];
    nWinTr = new uint[P.alignWindowsPerReadNmax];
    trArray = new Transcript[P.alignTranscriptsPerReadNmax];
    trArrayPointer =  new Transcript*[P.alignTranscriptsPerReadNmax];
    for (uint ii=0;ii<P.alignTranscriptsPerReadNmax;ii++)
        trArrayPointer[ii]= &(trArray[ii]);
    //read
    Read0 = new char*[2];
    Read0[0]  = new char [DEF_readSeqLengthMax+1];
    Read0[1]  = new char [DEF_readSeqLengthMax+1];
    Qual0 = new char*[2];
    Qual0[0]  = new char [DEF_readSeqLengthMax+1];
    Qual0[1]  = new char [DEF_readSeqLengthMax+1];
    readNameMates=new char* [P.readNmates];
    for (uint ii=0; ii<P.readNmates; ii++) {
        readNameMates[ii]=new char [DEF_readNameLengthMax];
    };
    readNameExtra.resize(P.readNmates);
    readName = readNameMates[0];
    Read1 = new char*[3];
    Read1[0]=new char[DEF_readSeqLengthMax+1]; Read1[1]=new char[DEF_readSeqLengthMax+1]; Read1[2]=new char[DEF_readSeqLengthMax+1];
    Qual1=new char*[2]; //modified QSs for scoring
    Qual1[0]=new char[DEF_readSeqLengthMax+1]; Qual1[1]=new char[DEF_readSeqLengthMax+1];
    //outBAM
    outBAMoneAlignNbytes = new uint [P.readNmates+2]; //extra piece for chimeric reads
    outBAMoneAlign = new char* [P.readNmates+2]; //extra piece for chimeric reads
    for (uint ii=0; ii<P.readNmates+2; ii++) {
        outBAMoneAlign[ii]=new char [BAMoutput_oneAlignMaxBytes];
    };
    resetN();
    //chim
    chunkOutChimJunction = new fstream;
};

void ReadAlign::resetN () {//reset resets the counters to 0 for a new read
    mapMarker=0;
    nA=0;nP=0;
    WC.clear();
    nTr=0;nTrMate=0;
    nUM[0]=0;nUM[1]=0;
    storedLmin=0; uniqLmax=0; uniqLmaxInd=0; multLmax=0; multLmaxN=0; multNminL=0; multNmin=0; multNmax=0; multNmaxL=0;
    chimN=0;

    for (uint ii=0; ii<P.readNmates; ii++) {
        maxScoreMate[ii]=0;
    };
};

