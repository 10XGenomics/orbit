#ifndef MMAPARRAY_DEF
#define MMAPARRAY_DEF

#include "IncludeDefine.h"

class MmapArray {
    private:
        int fd;

    public:
        uint file_mmap_length;
        char *file_mmap;
        char *prefix_mmap;
        char *last_page_mmap;


        MmapArray();
        ~MmapArray();
        int makeMmap(string filename, uint length, uint suffix_padding);
};

#endif
