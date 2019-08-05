#include "Genome.h"
#include "Parameters.h"
#include "ReadAlign.h"
#include "Transcriptome.h"
#include "Variation.h"
#include "orbit.h"

struct Aligner
{
    public:
        // p represents the command line parameters and everything
        // that is calculated firectly from them
        Parameters *p;

        // ra represents the ReadAlign object that is used to make any kind of
        // alignment queries
        ReadAlign *ra;

        // g is the information contained in the built genome index
        Genome *g;

        // isOriginal is true iff the aligner is initialized with init()
        // instead of init_clone(), and is used for deciding which members
        // originated in this instance and can be safely freed upon destruction
        int isOriginal;

        Aligner(int argInN, char* argIn[])
        {
            isOriginal = 1;
            p = new Parameters();
            p->inputParameters(argInN, argIn);
            g = new Genome(*p);
            g->genomeLoad();
            Transcriptome *mainTranscriptome = nullptr;
            g->Var = new Variation(*p, g->chrStart, g->chrNameIndex);
            ra = new ReadAlign(*p, *g, mainTranscriptome, 0);
        }

        // This constructor is used to construct clones of an existing Aligner
        // This allows multi-threaded alignment without each thread
        // constructing its own genome object
        Aligner(Aligner* og)
        {
            isOriginal = 0;
            p = og->p;
            g = og->g;
            Transcriptome *mainTranscriptome = nullptr;
            ra = new ReadAlign(*p, *g, mainTranscriptome, 0);
        }

        ~Aligner()
        {
            delete ra;
            if(isOriginal)
            {
                delete g;
                delete p;
            }
        }
};


const char* align_read(Aligner* a, char *Read1, char *Qual1, unsigned long long read_length)
{
    a->p->iReadAll++;
    a->ra->iRead++;
    a->p->readNmates = 1;
    a->ra->readNmates = 1;
    a->ra->Read0 = &Read1;
    a->ra->Qual0 = &Qual1;
    a->ra->readName = (char*)malloc(2);
    a->ra->readName[0] = 'a';
    a->ra->readName[1] = '\0';
    a->ra->readLength[0] = read_length;
    a->ra->readLengthOriginal[0] = read_length;
    
    int readStatus = a->ra->oneRead();
    a->ra->readName[1] = '\0';
    if(readStatus != 0)
    {
        return "";
    }
    const char* str = a->ra->outputAlignments();
    return str;
}

const char* align_read_pair(Aligner* a, char *Read1, char *Qual1, char *Read2, char *Qual2, unsigned long long read_length)
{
    a->p->iReadAll++;
    a->ra->iRead++;
    a->p->readNmates = 2;
    a->ra->readNmates = 2;
    a->ra->Read0 = &Read1;
    a->ra->Qual0 = &Qual1;
    strcpy(a->ra->Read0[1], Read2);
    strcpy(a->ra->Qual0[1], Qual2);
    a->ra->readName = (char*)malloc(2);
    a->ra->readName[0] = 'a';
    a->ra->readName[1] = '\0';
    a->ra->readLength[0] = read_length;
    a->ra->readLengthOriginal[0] = read_length;
    
    int readStatus = a->ra->oneRead();
    a->ra->readName[1] = '\0';
    if(readStatus != 0)
    {
        return "";
    }
    const char* str = a->ra->outputAlignments();
    return str;
}

Aligner* init_aligner_clone(Aligner* al)
{
    return new Aligner(al);
}

Aligner* init_aligner(int argc, char* argv[])
{
    return new Aligner(argc, argv);
}

void destroy_aligner(Aligner *a)
{
    delete a;
}

