#include <memory>  // for make_unique

#include "Genome.h"
#include "Parameters.h"
#include "ReadAlign.h"

#include "orbit.h"

using std::make_unique;
using std::unique_ptr;

struct StarRef {
    public:
        const unique_ptr<Parameters> p;
        const unique_ptr<Genome> g;
        StarRef(int argInN, const char* const argIn[]);
};

namespace {

unique_ptr<Parameters> make_parameters(int argInN, const char* const argIn[]) {
    unique_ptr<Parameters> pMut = make_unique<Parameters>();
    pMut->inputParameters(argInN, argIn);
    //pMut->readNmates = 1;
    return pMut;
}

unique_ptr<Genome> load_genome(Parameters& p) {
    unique_ptr<Genome> gMut = make_unique<Genome>(p);
    gMut->genomeLoad();
    gMut->Var = nullptr; //new Variation(*pMut, gMut->chrStart, gMut->chrNameIndex);
    return gMut;
}

unique_ptr<ReadAlign> make_ra(const StarRef *ref) {
    return make_unique<ReadAlign>(*(ref->p), *(ref->g), nullptr, 0);
}

}  // namespace

StarRef::StarRef(int argInN, const char* const argIn[])
    : p(make_parameters(argInN, argIn)), g(load_genome(*p))
{ }

struct Aligner final {
    private:
        const unique_ptr<StarRef> owned_ref;

    public:

        // A pointer to the built index containing parameters and the
        // reference genome
        const StarRef *ref;

        // ra represents the ReadAlign object that is used to make any kind of
        // alignment queries
        unique_ptr<ReadAlign> ra;

        explicit Aligner(const StarRef* r)
            : ref(r),
              ra(make_ra(ref))
        { }

        Aligner(int argInN, const char* const argIn[])
            : owned_ref(make_unique<StarRef>(argInN, argIn)),
              ref(owned_ref.get()),
              ra(make_ra(ref))
        { }

        // This constructor is used to construct clones of an existing Aligner
        // This allows multi-threaded alignment without each thread
        // constructing its own genome object
        explicit Aligner(const Aligner* og)
            : ref(og->ref),
              ra(make_ra(ref))
        { }
};


const char* align_read(Aligner* a, const char* read1Fastq) {
    a->ra->iRead++;
    a->ra->readNmates = 1;
    a->ra->readFastq[0] = read1Fastq;
    a->ra->readName = "a";
    int readStatus = a->ra->oneRead();
    if(readStatus != 0) {
        return nullptr;
    }
    const char* str = a->ra->outputAlignments();
    return str;
}

const char* align_read_pair(Aligner* a, const char* read1Fastq, const char* read2Fastq) {
    a->ra->iRead++;
    a->ra->readNmates = 2;
    a->ra->readFastq[0] = read1Fastq;
    a->ra->readFastq[1] = read2Fastq;
    a->ra->readName = "a";
    
    int readStatus = a->ra->oneRead();
    if(readStatus != 0) {
        return nullptr;
    }
    const char* str = a->ra->outputAlignments();
    return str;
}

Aligner* init_aligner_clone(const Aligner* al) {
    return new Aligner(al);
}

Aligner* init_aligner(int argc, const char* const argv[]) {
    return new Aligner(argc, argv);
}

const StarRef* init_star_ref(int argc, const char* const argv[]) {
    return new StarRef(argc, argv);
}

Aligner* init_aligner_from_ref(const StarRef* sr) {
    return new Aligner(sr);
}

void destroy_aligner(Aligner *a) {
    delete a;
}

void destroy_ref(const StarRef* sr) {
    delete sr;
}

