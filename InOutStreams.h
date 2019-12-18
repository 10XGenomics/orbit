#ifndef INOUTSTREAMS_DEF
#define INOUTSTREAMS_DEF

#include <ios>
#include <ostream>
#include "IncludeDefine.h"
#include SAMTOOLS_BGZF_H

template<class charT, class traits = std::char_traits<charT> >
class NullBuf: public std::basic_streambuf<charT, traits> {
    inline typename traits::int_type overflow(typename traits::int_type c) {
        return traits::not_eof(c);
    }
};

template<class charT = char, class traits = std::char_traits<charT> >
class NullStream: public std::basic_ostream<charT, traits> {
    public:
        inline NullStream():
            std::basic_ios<charT, traits>(&nullbuf),
            std::basic_ostream<charT, traits>(&nullbuf)
        { std::basic_ios<charT, traits>::init(&nullbuf); }
        inline void open(const char*, std::ios_base::openmode = std::ios_base::out) {}
        inline void close() {}

    private:
        NullBuf<charT, traits> nullbuf;
};

class InOutStreams {
    public:
    ostream *logStdOut, *outSAM;
    ofstream outSAMfile;
    BGZF *outBAMfileUnsorted, *outBAMfileCoord, *outQuantBAMfile;

    ofstream outChimSAM, outChimJunction, logFinal, outUnmappedReadsStream[MAX_N_MATES];
    ifstream readIn[MAX_N_MATES];

    //send logs to nothing
    NullStream<> logStdOutFile, logMain, logProgress;

    //compilation-optional streams
    ofstream outLocalChains;

    InOutStreams();
    ~InOutStreams();
};

#endif
