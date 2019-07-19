#ifdef __cplusplus
extern "C" {
#endif

    struct Aligner;
    const char* align_read(struct Aligner*, char*, char*, unsigned long long);
    struct Aligner* init_aligner(int, char*[]);
    void destroy_aligner(struct Aligner*);

#ifdef __cplusplus
}
#endif

