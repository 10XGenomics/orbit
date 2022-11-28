#include "IncludeDefine.h"
#include "Parameters.h"
#include "Transcript.h"
#include "ReadAlign.h"

ReadAlign::ReadAlign (const Parameters& Pin, const Genome &genomeIn, Transcriptome *TrIn, int iChunk)
                    : mapGen(genomeIn), readFastq{nullptr, nullptr}, P(Pin), chunkTr(TrIn)
{
    iRead = 0;
    readFilesIndex = 0;
    readNmates=P.readNmates;
    winBin[0] = make_unique<uintWinBin[]>(P.winBinN);
    winBin[1] = make_unique<uintWinBin[]>(P.winBinN);
    memset(winBin[0].get(),255,sizeof(winBin[0][0])*P.winBinN);
    memset(winBin[1].get(),255,sizeof(winBin[0][0])*P.winBinN);
    //RNGs
    rngMultOrder.seed(P.runRNGseed*(iChunk+1));
    rngUniformReal0to1=std::uniform_real_distribution<double> (0.0, 1.0);
    //transcriptome
    if ( P.quant.trSAM.yes ) {
        alignTrAll=make_unique<Transcript[]>(P.alignTranscriptsPerReadNmax);
    };
    //split
    splitR[0].resize(P.maxNsplit); splitR[1].resize(P.maxNsplit); splitR[2].resize(P.maxNsplit);
    //alignments
    PC=make_unique<uiPC[]>(P.seedPerReadNmax);
    WC=make_unique<uiWC[]>(P.alignWindowsPerReadNmax);
    nWA=make_unique<uint[]>(P.alignWindowsPerReadNmax);
    nWAP=make_unique<uint[]>(P.alignWindowsPerReadNmax);
    WALrec=make_unique<uint[]>(P.alignWindowsPerReadNmax);
    WlastAnchor=make_unique<uint[]>(P.alignWindowsPerReadNmax);

#ifdef COMPILE_FOR_LONG_READS
    swWinCov = new uint[P.alignWindowsPerReadNmax];
    scoreSeedToSeed = new intScore [P.seedPerWindowNmax*(P.seedPerWindowNmax+1)/2];
    scoreSeedBest = new intScore [P.seedPerWindowNmax];
    scoreSeedBestInd = new uint [P.seedPerWindowNmax];
    scoreSeedBestMM = new uint [P.seedPerWindowNmax];
    seedChain = new uint [P.seedPerWindowNmax];
#endif

    WA=make_unique<unique_ptr<uiWA[]>[]>(P.alignWindowsPerReadNmax);
    for (uint ii=0;ii<P.alignWindowsPerReadNmax;ii++)
        WA[ii]=make_unique<uiWA[]>(P.seedPerWindowNmax);
    WAincl = make_unique<bool[]>(P.seedPerWindowNmax);
    trAll = make_unique<Transcript**[]>(P.alignWindowsPerReadNmax+1);
    nWinTr = make_unique<uint[]>(P.alignWindowsPerReadNmax);
    trArray = make_unique<Transcript[]>(P.alignTranscriptsPerReadNmax);
    trArrayPointer = make_unique<Transcript*[]>(P.alignTranscriptsPerReadNmax);
    for (uint ii=0;ii<P.alignTranscriptsPerReadNmax;ii++)
        trArrayPointer[ii]= &(trArray[ii]);
    trInit = make_unique<Transcript>();
    //read
    readData = make_unique<std::array<std::array<char, DEF_readSeqLengthMax+1>, 9>>();
    Read0[0] = std::get<0>(*readData).data();
    Read0[1] = std::get<1>(*readData).data();
    Qual0[0] = std::get<2>(*readData).data();
    Qual0[1] = std::get<3>(*readData).data();
    readNameMates=make_unique<char*[]>(P.readNmates);
    readNameMatesData=make_unique<std::array<char, DEF_readNameLengthMax>[]>(P.readNmates);
    for (uint ii=0; ii<P.readNmates; ii++) {
        readNameMates[ii]=readNameMatesData[ii].data();
    };
    readNameExtra.resize(P.readNmates);
    readName = readNameMates[0];
    Read1[0]=std::get<4>(*readData).data(); Read1[1]=std::get<5>(*readData).data(); Read1[2]=std::get<6>(*readData).data();
    //modified QSs for scoring
    Qual1[0]=std::get<7>(*readData).data(); Qual1[1]=std::get<8>(*readData).data();
    //outBAM
    outBAMoneAlignNbytes = make_unique<uint[]>(P.readNmates+2); //extra piece for chimeric reads
    outBAMoneAlign = make_unique<char*[]>(P.readNmates+2); //extra piece for chimeric reads
    outBAMoneAlignData = make_unique<std::array<char, BAMoutput_oneAlignMaxBytes>[]>(P.readNmates+2); //extra piece for chimeric reads
    for (uint ii=0; ii<P.readNmates+2; ii++) {
        outBAMoneAlign[ii]=outBAMoneAlignData[ii].data();
    };
    resetN();
};

void ReadAlign::resetN () {//reset resets the counters to 0 for a new read
    mapMarker=0;
    nA=0;nP=0;nW=0;
    nTr=0;nTrMate=0;
    nUM[0]=0;nUM[1]=0;
    storedLmin=0; uniqLmax=0; uniqLmaxInd=0; multLmax=0; multLmaxN=0; multNminL=0; multNmin=0; multNmax=0; multNmaxL=0;
    chimN=0;

    for (uint ii=0; ii<P.readNmates; ii++) {
        maxScoreMate[ii]=0;
    };

//     for (uint ii=0;ii<P.alignTranscriptsPerReadNmax;ii++) trArrayPointer[ii]= &(trArray[ii]);

};

