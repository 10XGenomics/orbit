#ifndef MMAPARRAY_DEF
#define MMAPARRAY_DEF

#include "IncludeDefine.h"

class MmapArray {

    public:
        // Start of the mapped region
        char *mmap_addr;
        // Total length of the array
        size_t mmap_length;

        // Position where the first byte of the file is mapped
        char *file_mmap_addr;
        // Length of the array after the start of the file map
        size_t file_mmap_length;

        MmapArray();
        ~MmapArray();
        int makeMmap(string filename, size_t length, size_t suffix_padding);
};

#endif
