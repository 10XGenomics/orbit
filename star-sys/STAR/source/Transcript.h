#ifndef CODE_Transcript
#define CODE_Transcript

#include <array>

#include "IncludeDefine.h"
#include "Parameters.h"
#include "Variation.h"
#include "Genome.h"

class Transcript {
public:
    std::array<Exon, MAX_N_EXONS> exons; //coordinates of all exons: r-start, g-start, length
    std::array<std::array<uint, 2>, MAX_N_EXONS> shiftSJ; //shift of the SJ coordinates due to genomic micro-repeats
    std::array<int, MAX_N_EXONS> canonSJ; //canonicity of each junction
    std::array<uint8, MAX_N_EXONS> sjAnnot; //anotated or not
    std::array<uint8, MAX_N_EXONS> sjStr; //strand of the junction

    std::array<int, 3> intronMotifs;
    uint8 sjMotifStrand;

    uint nExons; //number of exons in the read transcript

    //variables from ReadAlign
    uint *readLengthOriginal, *readLength;
    uint Lread, readLengthPairOriginal;
    uint iRead; //read identifier
    uint readNmates;
    char *readName;

    int iFrag; //frag number of the transcript, if the transcript contains only one frag

    //loci
    uint rStart = 0; //read
    uint roStart = 0; //original read
    uint rLength = 0; //read length
    uint gStart = 0; //genomic start
    uint gLength = 0; //genomic length
    uint cStart; //chromosome start
    uint Chr,Str,roStr; //chromosome and strand and original read Strand

    bool primaryFlag = false;

    uint nMatch = 0;//min number of matches
    uint nMM = 0;//max number of mismatches
    uint mappedLength; //total mapped length, sum of lengths of all blocks(exons)

    uint extendL = 0; //extension length
    intScore maxScore = 0; //maximum Score

    uint nGap = 0;
    uint lGap = 0; //number of genomic gaps (>alignIntronMin) and their total length
    uint nDel = 0; //number of genomic deletions (ie genomic gaps)
    uint nIns = 0; //number of (ie read gaps)
    uint lDel = 0; //total genomic deletion length
    uint lIns = 0; //total genomic insertion length

    uint nUnique = 0; //number of anchor pieces in the alignment
    uint nAnchor = 0; //number of unique pieces in the alignment

    vector <int32> varInd;
    vector <int32> varGenCoord, varReadCoord ;
    vector <char> varAllele;

    void reset(); //reset to 0
    void add(const Transcript*) noexcept; // add
    intScore alignScore(const char **Read1, const char *G, const Parameters &P);
    int variationAdjust(const Genome &mapGen, const char *R);
    string generateCigarP() const; //generates CIGAR
    void peOverlapSEtoPE(const uint* mSta, const Transcript &t);
    bool extractSpliceJunctions(vector<array<uint64,2>> &sjOut) const;

    Transcript() = default;
};

#endif
