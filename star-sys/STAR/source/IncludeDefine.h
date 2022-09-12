#ifndef INCLUDEDEFINE_DEF
#define INCLUDEDEFINE_DEF

//standard libs
#include <algorithm>
#include <cstring>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <iomanip>
#include <vector>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <cerrno>
#include <limits>
#include <cstdint>

#include "VERSION"

#define ERROR_OUT string ( __FILE__ ) +":"+ to_string ( (uint) __LINE__ ) +":"+ string ( __FUNCTION__ )

//external libs
#define SAMTOOLS_BGZF_H "htslib/htslib/bgzf.h"
#define SAMTOOLS_SAM_H  "htslib/htslib/sam.h"

using namespace std;

#if defined(__mips__) && !defined(SHM_NORESERVE)
#define SHM_NORESERVE 010000
#endif

typedef int8_t int8;
typedef uint8_t uint8;

#define uint unsigned long long
#define sint signed long long
#define uint64 unsigned long long
#define uint32 unsigned int
#define uint16 unsigned short int
#define uchar unsigned char
#define int64 long long
#define int32 int

// this is gcc extension, may need to redefine for other compilers
#define uint128 __uint128_t

#define GENOME_spacingChar 5

#define uintWinBin unsigned short
#define uintWinBinMax numeric_limits<uint16>::max()


#define intSWscore int
#define intScore int

#define scoreMatch 1


//cleaned
//output
#define BAMoutput_oneAlignMaxBytes 100000


//SAM attributes
#define ATTR_NH 1
#define ATTR_HI 2
#define ATTR_AS 3
#define ATTR_NM 4
#define ATTR_MD 5
#define ATTR_nM 6
#define ATTR_jM 7
#define ATTR_jI 8
#define ATTR_XS 9
#define ATTR_RG 10
#define ATTR_vG 11
#define ATTR_vA 12
#define ATTR_vW 13
#define ATTR_ch 14
#define ATTR_MC 15
#define ATTR_rB 16
#define ATTR_CR 17
#define ATTR_CY 18
#define ATTR_UR 19
#define ATTR_UY 20

#define MAX_N_EXONS 20

//input reads
#define MAX_N_MATES 2
#define DEF_readNameLengthMax 50000
    #define DEF_readSeqLengthMax 650

#if (DEF_readNameLengthMax > DEF_readSeqLengthMax)
        #define DEF_readNameSeqLengthMax DEF_readNameLengthMax
#else
        #define DEF_readNameSeqLengthMax DEF_readSeqLengthMax
#endif

#define EXIT_CODE_BUG 101
#define EXIT_CODE_PARAMETER 102
#define EXIT_CODE_RUNTIME 103
#define EXIT_CODE_INPUT_FILES 104
#define EXIT_CODE_GENOME_FILES 105
#define EXIT_CODE_MEMORY_ALLOCATION 108
#define EXIT_CODE_FILE_OPEN 109
#define EXIT_CODE_FILE_WRITE 110

//cleaned-end


//exit codes
#define EXIT_createExtendWindowsWithAlign_TOO_MANY_WINDOWS 101

#define SJ_MOTIF_SIZE 7 //number of recorded SJ motifs
#define SJ_SAM_AnnotatedMotifShift 20

#define EXTEND_ORDER 1 //1-first extend to the 5' of the read, then 3'; 2- first extend to the left, then to the right

#define MARK_FRAG_SPACER_BASE 11
#define MAX_N_CHIMERAS 5
#define MAX_N_MULTMAP 100000 //max number of multiple mappers
#define MAX_SJ_REPEAT_SEARCH 255 //max length of a repeat to search around a SJ



struct Window {
    uint Str, Chr, gStart, gEnd;
};


struct Exon {
    uint R, G, L, iFrag, sjA;
};



#define MARKER_ALL_PIECES_EXCEED_seedMultimapNmax 999901 //marks the reads that map too many time, more than seedMultimapNmax
#define MARKER_NO_GOOD_WINDOW 999903 //did not find any good windows
#define MARKER_NO_GOOD_PIECES 999904
#define MARKER_TOO_MANY_ANCHORS_PER_WINDOW 999905
#define MARKER_READ_TOO_SHORT 999910

struct PC {
uint rStart, Length, Dir, Nrep, SAstart, SAend, iFrag;
};

struct uiWA {
    uint Length, rStart, gStart, Nrep, Anchor;
    // Which fragment is this read from? (typically 0 for single end data,
    // can be > for paired end.
    uint iFrag;
    // sjA is (uint) -1 if the alignment is not from a splice junction, or
    // else is the index of the splice junction sequence this represents
    // (first splice junction in suffix array is 0);
    uint sjA;

    public:
        uiWA() = default;
        uiWA(uint Length, uint rStart, uint gStart, uint Nrep, uint Anchor, uint iFrag, uint sjA) :
        Length(Length),
        rStart(rStart),
        gStart(gStart),
        Nrep(Nrep),
        Anchor(Anchor),
        iFrag(iFrag),
        sjA(sjA)
        {}
};

// debugging
#if defined DEBUG
    #define DEBUG_stitch
    #define DEBUG_Nread 200000
    #define DEBUG_NreadStart 1
    #define DEBUG_extend
#endif

#endif
