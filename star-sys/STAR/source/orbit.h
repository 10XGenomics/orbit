/*
 * Orbit is a wrapper around the key functionality of STAR.  It supports
 * querying a previously built genome index with single reads and read pairs
 * and obtaining BAM records with their alignments to the index.
 */

#ifdef __cplusplus
extern "C" {
#endif
    // StarRef: A separate struct for the parameters and genome information.
    // This is so that this can be built in the main thread and declared as
    // const to ensure safe access by the other threads
    struct StarRef;

    // Aligner: all constructed alignment object which is used to align
    // individual reads/read pairs through the functions below
    struct Aligner;

    // align_read: align an individual read and get a string of BAM records
    const char* align_read(struct Aligner*, const char*);
    
    // align_read_pair: align a pair of reads and get a string of BAM records
    const char* align_read_pair(struct Aligner*, const char*, const char*);

    // init_aligner_clone: create an aligner from the same reference as an
    // existing aligner, sharing key structures with it and saving memory in
    // multi-threaded applications
    // NOTE: this function is deprecated in favor of a StarRef-based workflow
    struct Aligner* init_aligner_clone(const struct Aligner*);

    // init_aligner: initialize an aligner given the array of parameters which
    // would be passed to STAR
    struct Aligner* init_aligner(int, const char* const[]);

    // init_star_ref: build a star reference with a given set of arguments
    const struct StarRef* init_star_ref(int, const char* const[]);

    // init_aligner_from_ref takes a StarRef struct with an already built
    // genome and builds an aligner around it
    struct Aligner* init_aligner_from_ref(const struct StarRef*);

    // destroy_aligner: frees the memory occupied by an aligner
    void destroy_aligner(struct Aligner*);

    //destroy_ref: frees the memory occupied by a reference
    void destroy_ref(const struct StarRef*);

#ifdef __cplusplus
}
#endif

