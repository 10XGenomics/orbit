#ifndef CODE_ChimericSegment
#define CODE_ChimericSegment

#include "IncludeDefine.h"
#include "Parameters.h"
#include "Transcript.h"
#include "ParametersChimeric.h"

class ChimericSegment
{//
    public:
        const Parameters &P;
        const ParametersChimeric &pCh;

        Transcript  &align;     //alignment
        uint roS,roE,str;     //start/end/strand in original read coordinates

        ChimericSegment(const Parameters &Pin, Transcript &alignIn); //allocate
        bool segmentCheck();//check if chimeric segment is good
    private:
};

#endif
