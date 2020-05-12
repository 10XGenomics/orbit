#include "Genome.h"
#include "Parameters.h"
#include "ReadAlign.h"
#include "Transcriptome.h"

#include "orbit.h"

struct StarRef
{
    public:
        const Parameters* p;
        const Genome* g;
        StarRef(int argInN, char* argIn[])
        {
            Parameters* pMut = new Parameters();
            pMut->inputParameters(argInN, argIn);
            //pMut->readNmates = 1;
            Genome* gMut = new Genome(*pMut);
            gMut->genomeLoad();
            gMut->Var = NULL; //new Variation(*pMut, gMut->chrStart, gMut->chrNameIndex);
            p = pMut;
            g = gMut;
        }
        ~StarRef()
        {
            delete g;
            delete p;
        }
};

struct Aligner
{
    public:

        // A pointer to the built index containing parameters and the
        // reference genome
        const StarRef *ref;

        // ra represents the ReadAlign object that is used to make any kind of
        // alignment queries
        ReadAlign *ra;

        // isOriginal is true iff the aligner is initialized with init()
        // instead of init_clone(), and is used for deciding which members
        // originated in this instance and can be safely freed upon destruction
        int isOriginal;

        Aligner(const StarRef* r)
        {
            isOriginal = 0;
            ref = r;
            Transcriptome *mainTranscriptome = nullptr;
            ra = new ReadAlign(*(ref->p), *(ref->g), mainTranscriptome, 0);
        }

        Aligner(int argInN, char* argIn[])
        {
            isOriginal = 1;
            ref = new StarRef(argInN, argIn);
            Transcriptome *mainTranscriptome = nullptr;
            ra = new ReadAlign(*(ref->p), *(ref->g), mainTranscriptome, 0);
        }

        // This constructor is used to construct clones of an existing Aligner
        // This allows multi-threaded alignment without each thread
        // constructing its own genome object
        Aligner(const Aligner* og)
        {
            isOriginal = 0;
            ref = og->ref;
            Transcriptome *mainTranscriptome = nullptr;
            ra = new ReadAlign(*(ref->p), *(ref->g), mainTranscriptome, 0);
        }

        ~Aligner()
        {
            delete ra;
            if(isOriginal)
            {
                delete ref;
            }
        }
};


const char* align_read(Aligner* a, const char* read1Fastq)
{
    static char qname[] = "a";
    a->ra->iRead++;
    a->ra->readNmates = 1;
    a->ra->readFastq[0] = read1Fastq;
    a->ra->readName = qname;
    int readStatus = a->ra->oneRead();
    a->ra->readName[1] = '\0';
    if(readStatus != 0)
    {
        return nullptr;
    }
    const char* str = a->ra->outputAlignments();
    return str;
}

const char* align_read_pair(Aligner* a, const char* read1Fastq, const char* read2Fastq)
{
    static char qname[] = "a";
    a->ra->iRead++;
    a->ra->readNmates = 2;
    a->ra->readFastq[0] = read1Fastq;
    a->ra->readFastq[1] = read2Fastq;
    a->ra->readName = qname;
    
    int readStatus = a->ra->oneRead();
    a->ra->readName[1] = '\0';
    if(readStatus != 0)
    {
        return nullptr;
    }
    const char* str = a->ra->outputAlignments();
    return str;
}

Aligner* init_aligner_clone(const Aligner* al)
{
    return new Aligner(al);
}

Aligner* init_aligner(int argc, char* argv[])
{
    return new Aligner(argc, argv);
}

const StarRef* init_star_ref(int argc, char* argv[])
{
    return new StarRef(argc, argv);
}

Aligner* init_aligner_from_ref(const StarRef* sr)
{
    return new Aligner(sr);
}

void destroy_aligner(Aligner *a)
{
    delete a;
}

void destroy_ref(const StarRef* sr)
{
    delete sr;
}

