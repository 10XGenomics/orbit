#include "Genome.h"
#include "Parameters.h"
#include "SuffixArrayFuns.h"
#include "PackedArray.h"
#include "ErrorWarning.h"
#include "streamFuns.h"
#include "genomeScanFastaFiles.h"

#include <time.h>
#include <cmath>
#include <unistd.h>
#include <sys/stat.h>


Genome::Genome (Parameters &Pin ): pGe(Pin.pGe), P(Pin) {


    sjdbOverhang = pGe.sjdbOverhang; //will be re-defined later if another value was used for the generated genome
    sjdbLength = pGe.sjdbOverhang==0 ? 0 : pGe.sjdbOverhang*2+1;
};

Genome::~Genome()
{
    freeMemory();
}

void Genome::freeMemory(){//free big chunks of memory used by genome and suffix array

    if (pGe.gLoad=="NoSharedMemory") {//can deallocate only for non-shared memory
        SA.deallocateArray();
        SApass2.deallocateArray();
        SAi.deallocateArray();
    };
};

uint Genome::OpenStream(string name, ifstream & stream, uint size)
{
    stream.open((pGe.gDir+ "/" +name).c_str(), ios::binary);
    if (!stream.good()) {
        ostringstream errOut;
        errOut << "EXITING because of FATAL ERROR: could not open genome file: "<< pGe.gDir << "/" << name <<"\n";
        errOut << "SOLUTION: check that the path to genome files, specified in --genomeDir is correct and the files are present, and have user read permissions\n" <<flush;
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_GENOME_FILES, P);
    };


    if (size>0) {
        P.inOut->logMain << name << ": size given as a parameter = " << size <<"\n";
    } else {
        P.inOut->logMain << "Checking " << name << " size";
        stream.seekg (0, ios::end);
        int64 size1 = stream.tellg();
        if (size1<=0) {
            ostringstream errOut;
            errOut << "EXITING because of FATAL ERROR: failed reading from genome file: "<< pGe.gDir << "/" << name <<"\n";
            errOut << "SOLUTION: re-generate the genome index\n";
            exitWithError(errOut.str(),std::cerr, P.inOut->logMain, 1, P);
        };
        size=(uint) size1;
        stream.clear();
        stream.seekg (0, ios::beg);
        P.inOut->logMain << "file size: "<< size <<" bytes; state: good=" <<stream.good()\
                <<" eof="<<stream.eof()<<" fail="<<stream.fail()<<" bad="<<stream.bad()<<"\n"<<flush;
    };

    return size;
};


void Genome::genomeLoad(){//allocate and load Genome

    time_t rawtime;
    time ( &rawtime );
    *(P.inOut->logStdOut) << timeMonthDayTime(rawtime) << " ..... loading genome\n" <<flush;

    uint L=200,K=6;

    Parameters P1;

    //some initializations before reading the parameters
    GstrandBit=0;

    ifstream parFile((pGe.gDir+("/genomeParameters.txt")).c_str());
    if (parFile.good()) {
        P.inOut->logMain << "Reading genome generation parameters:\n";

        //read genome internal parameters
        while (parFile.good()) {
            string word1;
            parFile >> word1;
            if (word1=="###") {
                parFile >> word1;
                if (word1=="GstrandBit") {
                    uint gsb1=0;
                    parFile >> gsb1;
                    GstrandBit=(uint8) gsb1;
                    P.inOut->logMain << "### GstrandBit=" << (uint) GstrandBit <<"\n";
                } else {
                    P.inOut->logMain << "### " <<word1;
                    getline(parFile,word1);
                    P.inOut->logMain <<word1<<"\n";
                };
            };
        };
        parFile.clear();
        parFile.seekg(0,ios::beg);//rewind

        // free P1.inOut first to avoid resource leak
        delete P1.inOut;
        P1.inOut = P.inOut;
        P1.scanAllLines(parFile,3,-1);
        parFile.close();
    } else {
        ostringstream errOut;
        errOut << "EXITING because of FATAL ERROR: could not open genome file "<< pGe.gDir+("/genomeParameters.txt") << endl;
        errOut << "SOLUTION: check that the path to genome files, specified in --genomeDir is correct and the files are present, and have user read permsissions\n" <<flush;
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_GENOME_FILES, P);
    };

    //check genome version
    if (P1.versionGenome.size()==0) {//
        ostringstream errOut;
        errOut << "EXITING because of FATAL ERROR: read no value for the versionGenome parameter from genomeParameters.txt file\n";
        errOut << "SOLUTION: please re-generate genome from scratch with the latest version of STAR\n";
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_GENOME_FILES, P);
    } else if (P1.versionGenome == P.versionGenome || P1.versionGenome == "20201") {//
        P.inOut->logMain << "Genome version is compatible with current STAR\n";
    } else {
        ostringstream errOut;
        errOut << "EXITING because of FATAL ERROR: Genome version: " << P1.versionGenome << " is INCOMPATIBLE with running STAR version: "<< STAR_VERSION <<"\n";
        errOut << "SOLUTION: please re-generate genome from scratch with running version of STAR, or with version: " << P.versionGenome <<"\n";
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_GENOME_FILES, P);
    };

    //find chr starts from files
    chrInfoLoad();

    //check if sjdbInfo.txt exists => genome was generated with junctions
    bool sjdbInfoExists=false;
    struct stat sjdb1;
    if ( stat( (pGe.gDir+"/sjdbInfo.txt").c_str(), &sjdb1) == 0 )
    {//file exists
        sjdbInfoExists=true;
    };

    if ( P.sjdbInsert.yes && sjdbInfoExists && P1.pGe.sjdbInsertSave=="")
    {//if sjdbInsert, and genome had junctions, and genome is old - it should be re-generated with new STAR
        ostringstream errOut;
        errOut << "EXITING because of FATAL ERROR: old Genome is INCOMPATIBLE with on the fly junction insertion\n";
        errOut << "SOLUTION: please re-generate genome from scratch with the latest version of STAR\n";
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_GENOME_FILES, P);
    };

    //record required genome parameters in P
    pGe.gSAindexNbases=P1.pGe.gSAindexNbases;
    pGe.gChrBinNbits=P1.pGe.gChrBinNbits;
    genomeChrBinNbases=1LLU<<pGe.gChrBinNbits;
    pGe.gSAsparseD=P1.pGe.gSAsparseD;

    if (P1.pGe.gFileSizes.size()>0)
    {//genomeFileSize was recorded in the genomeParameters file, copy the values to P
        pGe.gFileSizes = P1.pGe.gFileSizes;
    };

    if (P.parArray.at(pGe.sjdbOverhang_par)->inputLevel==0 && P1.pGe.sjdbOverhang>0)
    {//if --sjdbOverhang was not defined by user and it was defined >0 at the genome generation step, then use pGe.sjdbOverhang from the genome generation step
        pGe.sjdbOverhang=P1.pGe.sjdbOverhang;
        P.inOut->logMain << "--sjdbOverhang = " << pGe.sjdbOverhang << " taken from the generated genome\n";
    } else if (sjdbInfoExists && P.parArray.at(pGe.sjdbOverhang_par)->inputLevel>0 && pGe.sjdbOverhang!=P1.pGe.sjdbOverhang)
    {//if pGe.sjdbOverhang was defined at the genome generation step,the mapping step value has to agree with it
        ostringstream errOut;
        errOut << "EXITING because of fatal PARAMETERS error: present --sjdbOverhang="<<pGe.sjdbOverhang << " is not equal to the value at the genome generation step ="<< P1.pGe.sjdbOverhang << "\n";
        errOut << "SOLUTION: \n" <<flush;
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_GENOME_FILES, P);
    };

    sjdbOverhang = pGe.sjdbOverhang;
    sjdbLength = pGe.sjdbOverhang==0 ? 0 : pGe.sjdbOverhang*2+1;

    P.inOut->logMain << "Started loading the genome: " << asctime (localtime ( &rawtime ))<<"\n"<<flush;

    ifstream GenomeIn, SAin, SAiIn;

    if (pGe.gFileSizes.size() < 2)
    {//no size info available
        pGe.gFileSizes.push_back(0);
        pGe.gFileSizes.push_back(0);
    };
    nGenome = OpenStream("Genome",GenomeIn,pGe.gFileSizes.at(0));
    nSAbyte = OpenStream("SA",SAin,pGe.gFileSizes.at(1));
    OpenStream("SAindex",SAiIn,1); //we do not need SAiIn siz, using a dummy value here to prevent from reading its size from the disk

    uint SAiInBytes=0;
    SAiInBytes += fstreamReadBig(SAiIn,(char*) &pGe.gSAindexNbases, sizeof(pGe.gSAindexNbases));

    genomeSAindexStart = new uint[pGe.gSAindexNbases+1];

    // Note this data from SAiIn is small, don't bother mmaping it.
    SAiInBytes += fstreamReadBig(SAiIn,(char*) genomeSAindexStart, sizeof(genomeSAindexStart[0])*(pGe.gSAindexNbases+1));
    nSAi=genomeSAindexStart[pGe.gSAindexNbases];
    P.inOut->logMain << "Read from SAindex: pGe.gSAindexNbases=" << pGe.gSAindexNbases <<"  nSAi="<< nSAi <<endl;

    /////////////////////////////////// at this point all array sizes should be known: calculate packed array lengths
    if (GstrandBit==0) {//not defined before
        GstrandBit = (uint) floor(log(nGenome)/log(2))+1;
        if (GstrandBit<32) GstrandBit=32; //TODO: use simple access function for SA
    };


    GstrandMask = ~(1LLU<<GstrandBit);
    nSA=(nSAbyte*8)/(GstrandBit+1);
    SA.defineBits(GstrandBit+1,nSA);

    SAiMarkNbit=GstrandBit+1;
    SAiMarkAbsentBit=GstrandBit+2;

    SAiMarkNmaskC=1LLU << SAiMarkNbit;
    SAiMarkNmask=~SAiMarkNmaskC;
    SAiMarkAbsentMaskC=1LLU << SAiMarkAbsentBit;
    SAiMarkAbsentMask=~SAiMarkAbsentMaskC;

    SAi.defineBits(GstrandBit+3,nSAi);

    P.inOut->logMain << "nGenome=" << nGenome << ";  nSAbyte=" << nSAbyte <<endl<< flush;
    P.inOut->logMain <<"GstrandBit="<<int(GstrandBit)<<"   SA number of indices="<<nSA<<endl<<flush;


    genomeInsertL=0;
    genomeInsertChrIndFirst=nChrReal;

    //load genome
    P.inOut->logMain <<"Genome file size: "<<nGenome <<" bytes; state: good=" <<GenomeIn.good()\
            <<" eof="<<GenomeIn.eof()<<" fail="<<GenomeIn.fail()<<" bad="<<GenomeIn.bad()<<"\n"<<flush;
    P.inOut->logMain <<"Loading Genome ... " << flush;

    int res = mmapGenome.initMmap((pGe.gDir+ "/" + "Genome"), nGenome, L);
    if (res != 0) {
        ostringstream errOut;
        errOut <<"EXITING: got error in mmap: " << res;
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_MEMORY_ALLOCATION, P);
    }
    // G points to the start of genome
    G = mmapGenome.file_mmap_addr;

    // G1 points to the L bytes of padding before the genome. 
    // There is always one page of padding available.
    G1 = mmapGenome.file_mmap_addr - L;

    GenomeIn.close();

    for (uint ii=0;ii<L;ii++) {// attach a tail with the largest symbol
        G1[ii]=K-1;
        G[nGenome+ii]=K-1;
    };

    //load SAs
    P.inOut->logMain <<"SA file size: "<<SA.lengthByte <<" bytes; state: good=" <<SAin.good()\
            <<" eof="<<SAin.eof()<<" fail="<<SAin.fail()<<" bad="<<SAin.bad()<<"\n"<<flush;
    P.inOut->logMain <<"Loading SA ... " << flush;

    // Load suffix array with mmap
    res = mmapSA.initMmap((pGe.gDir+ "/" + "SA"), SA.lengthByte, 0);
    if (res != 0) {
        ostringstream errOut;
        errOut <<"EXITING: got error in mmap: " << res;
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_MEMORY_ALLOCATION, P);
    }

    SA.pointArray(mmapSA.file_mmap_addr);
    SAin.close();

    P.inOut->logMain <<"Loading SAindex ... " << flush;

    size_t SAindexOffset = SAiIn.tellg();
    res = mmapSAi.initMmap((pGe.gDir+ "/" + "SAindex"), SAindexOffset + SAi.lengthByte, 0);
    if (res != 0) {
        ostringstream errOut;
        errOut <<"EXITING: got error in mmap: " << res;
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_MEMORY_ALLOCATION, P);
    }

    // SAi starts past the begining of the file.
    SAi.pointArray(mmapSAi.file_mmap_addr + SAindexOffset);

    SAiIn.close();

    time ( &rawtime );
    P.inOut->logMain << "Finished ljk loading the genome: " << asctime (localtime ( &rawtime )) <<"\n"<<flush;

    insertSequences();

    chrBinFill();

    //splice junctions database
    if (nGenome==chrStart[nChrReal]) {//no sjdb
        sjdbN=0;
        sjGstart=chrStart[nChrReal]+1; //not sure why I need that
    } else {//there are sjdb chromosomes
        ifstream sjdbInfo((pGe.gDir+"/sjdbInfo.txt").c_str());
        if (sjdbInfo.fail()) {
            ostringstream errOut;
            errOut << "EXITING because of FATAL error, could not open file " << (pGe.gDir+"/sjdbInfo.txt") <<"\n";
            errOut << "SOLUTION: check that the path to genome files, specified in --genomeDir is correct and the files are present, and have user read permsissions\n" <<flush;
            exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        };


        sjdbInfo >> sjdbN >> pGe.sjdbOverhang;
        P.inOut->logMain << "Processing splice junctions database sjdbN=" <<sjdbN<<",   pGe.sjdbOverhang=" <<pGe.sjdbOverhang <<" \n";

        sjChrStart=nChrReal;
        sjGstart=chrStart[sjChrStart];

        //fill the sj-db to genome translation array
        sjDstart=new uint [sjdbN];
        sjAstart=new uint [sjdbN];
        sjdbStart=new uint [sjdbN];
        sjdbEnd=new uint [sjdbN];

        sjdbMotif=new uint8 [sjdbN];
        sjdbShiftLeft=new uint8 [sjdbN];
        sjdbShiftRight=new uint8 [sjdbN];
        sjdbStrand=new uint8 [sjdbN];

        for (uint ii=0;ii<sjdbN;ii++) {//get the info about junctions from sjdbInfo.txt
            {
                uint16 d1,d2,d3,d4;
                sjdbInfo >> sjdbStart[ii] >> sjdbEnd[ii] >> d1 >> d2 >> d3 >> d4;
                sjdbMotif[ii]      = (uint8) d1;
                sjdbShiftLeft[ii]  = (uint8) d2;
                sjdbShiftRight[ii] = (uint8) d3;
                sjdbStrand[ii] = (uint8) d4;
            };
            sjDstart[ii]   = sjdbStart[ii]  - pGe.sjdbOverhang;
            sjAstart[ii]   = sjdbEnd[ii] + 1;
            if (sjdbMotif[ii]==0) {//shinon-canonical junctions back to their true coordinates
                sjDstart[ii] += sjdbShiftLeft[ii];
                sjAstart[ii] += sjdbShiftLeft[ii];
            };
        };
    };

    //check and redefine some parameters
    //max intron size
    if (P.alignIntronMax==0 && P.alignMatesGapMax==0) {
        P.inOut->logMain << "alignIntronMax=alignMatesGapMax=0, the max intron size will be approximately determined by (2^winBinNbits)*winAnchorDistNbins=" \
                << (1LLU<<P.winBinNbits)*P.winAnchorDistNbins <<endl;
    } else {
        //redefine winBinNbits
        P.winBinNbits = (uint) floor( log2( max( max(4LLU,P.alignIntronMax), (P.alignMatesGapMax==0 ? 1000LLU : P.alignMatesGapMax) ) /4 ) + 0.5);
        P.winBinNbits = max( P.winBinNbits, (uint) floor(log2(nGenome/40000+1)+0.5) );
        //ISSUE - to be fixed in STAR3: if alignIntronMax>0 but alignMatesGapMax==0, winBinNbits will be defined by alignIntronMax
        P.inOut->logMain << "To accommodate alignIntronMax="<<P.alignIntronMax<<" redefined winBinNbits="<< P.winBinNbits <<endl;

    };

    if (P.winBinNbits > pGe.gChrBinNbits) {
       P.inOut->logMain << "winBinNbits=" <<P.winBinNbits <<" > " << "pGe.gChrBinNbits=" << pGe.gChrBinNbits << "   redefining:\n";
       P.winBinNbits=pGe.gChrBinNbits;
       P.inOut->logMain << "winBinNbits=" <<P.winBinNbits <<endl;
    };


    if (P.alignIntronMax==0 && P.alignMatesGapMax==0) {
    } else {
        //redefine winFlankNbins,winAnchorDistNbins
        P.winFlankNbins=max(P.alignIntronMax,P.alignMatesGapMax)/(1LLU<<P.winBinNbits)+1;
        P.winAnchorDistNbins=2*P.winFlankNbins;
        P.inOut->logMain << "To accommodate alignIntronMax="<<P.alignIntronMax<<" and alignMatesGapMax="<<P.alignMatesGapMax<<\
                ", redefined winFlankNbins="<<P.winFlankNbins<<" and winAnchorDistNbins="<<P.winAnchorDistNbins<<endl;
    };

    P.winBinChrNbits=pGe.gChrBinNbits-P.winBinNbits;
    P.winBinN = nGenome/(1LLU << P.winBinNbits)+1;//this may be chenaged later

    // reset P1.inOut here to avoid double-free (of P.inOut later)
    P1.inOut = nullptr;
};

//////////////////////////////////////////////////////////////////////////////////////////
void Genome::chrInfoLoad() {//find chrStart,Length,nChr from Genome G

    //load chr names
    ifstream chrStreamIn ( (pGe.gDir+"/chrName.txt").c_str() );
    if (chrStreamIn.fail()) {
        ostringstream errOut;
        errOut << "EXITING because of FATAL error, could not open file " << (pGe.gDir+"/chrName.txt") <<"\n";
        errOut << "SOLUTION: re-generate genome files with STAR --runMode genomeGenerate\n";
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
    };

    char chrInChar[1000];

    while (chrStreamIn.good()) {
        string chrIn;
        chrStreamIn.getline(chrInChar,1000);
        chrIn=chrInChar;
        if (chrIn=="") break;
        chrName.push_back(chrIn);
    };
    chrStreamIn.close();
    nChrReal=chrName.size();

    P.inOut->logMain << "Number of real (reference) chromosomes= " << nChrReal <<"\n"<<flush;
    chrStart.resize(nChrReal+1);
    chrLength.resize(nChrReal);

    //load chr lengths
    chrStreamIn.open( (pGe.gDir+"/chrLength.txt").c_str() );
    if (chrStreamIn.fail()) {
        ostringstream errOut;
        errOut << "EXITING because of FATAL error, could not open file " << (pGe.gDir+"/chrLength.txt") <<"\n";
        errOut << "SOLUTION: re-generate genome files with STAR --runMode genomeGenerate\n";
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
    };

    for  (uint ii=0;ii<nChrReal;ii++) {
        chrStreamIn >> chrLength[ii];
    };
    chrStreamIn.close();

    //load chr starts
    chrStreamIn.open( (pGe.gDir+"/chrStart.txt").c_str() );
    if (chrStreamIn.fail()) {
        ostringstream errOut;
        errOut << "EXITING because of FATAL error, could not open file " << (pGe.gDir+"/chrStart.txt") <<"\n";
        errOut << "SOLUTION: re-generate genome files with STAR --runMode genomeGenerate\n";
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
    };

    for  (uint ii=0;ii<=nChrReal;ii++) {
        chrStreamIn >> chrStart[ii];
    };
    chrStreamIn.close();

    //log
    for (uint ii=0; ii<nChrReal;ii++) {
        P.inOut->logMain << ii+1 <<"\t"<< chrName[ii] <<"\t"<<chrLength[ii]<<"\t"<<chrStart[ii]<<"\n"<<flush;
        chrNameIndex[chrName[ii]]=ii;
    };
};

//////////////////////////////////////////////////////////
void Genome::chrBinFill() {
    chrBinN = chrStart[nChrReal]/genomeChrBinNbases+1;
    chrBin = new uint [chrBinN];
    for (uint ii=0, ichr=1; ii<chrBinN; ++ii) {
        if (ii*genomeChrBinNbases>=chrStart[ichr]) ichr++;
        chrBin[ii]=ichr-1;
    };
};
