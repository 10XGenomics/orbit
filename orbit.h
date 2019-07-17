#include <sstream>
#include <string>

#include "Genome.h"
#include "Parameters.h"
#include "ReadAlign.h"
#include "Transcriptome.h"
#include "Variation.h"

extern "C" {
    struct Aligner;
    const char* align_read(Aligner*, char*, char*, unsigned long long);
    Aligner* init_aligner(int, char*[]);
    void destroy_aligner(Aligner*);
}
