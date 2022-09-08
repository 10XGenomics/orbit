#include <iostream>
#include "orbit.h"

int main() {
    char* input[16] = {
    "STAR",
    "--genomeDir",
    "/Users/nigel.delaney/refdata/refdata-gex-GRCh38-2020-A/star",
    "--runThreadN",
    "1",
    "--readNameSeparator",
    "space",
    "--outSAMunmapped",
    "Within",
    "KeepPairs",
    "--outSAMtype",
    "SAM",
    "--outStd",
    "SAM",
    "--outSAMorder",
    "PairedKeepInputOrder",
    };
    const struct StarRef* ref = init_star_ref(16, input);
    Aligner* al = init_aligner_from_ref(ref);
    const char* read = "@a\nGGATCACTTGAGGTCAGGAGTCCAAGACCAGCCTGGCCAACATGGTGAAACCCCCATCTCTACTCAAAATCCAAAAATTAGCCAGGCATA\n+\nJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ\n";
    //const char* read = "@a\nGGATCACTTGAGGTCAGGAGTCCAAGACCAGCCTGGCCAACATGGTGAAACCCNCATCTCTACTCAAAATCCAAAAATTAGCCAGGCATA\n+\nJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ\n";
    for(int i =0; i< 100000; i++) {
        align_read(al, read);
    }
    return 0;
}
