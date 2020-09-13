#ifndef MMAPARRAY_DEF
#define MMAPARRAY_DEF

#include "IncludeDefine.h"

class MmapArray {
    private:
        int fd;

    public:
        uint file_mmap_length;
        void *file_mmap;
        void *prefix_mmap;


        MmapArray();
        ~MmapArray();
        int makeMmap(string filename, uint length, uint suffix_padding);
};

#endif
