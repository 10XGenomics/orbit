#ifndef MMAPARRAY_H_
#define MMAPARRAY_H_

#include <string>
#include <cstddef>

class MmapArray {
 private:
   MmapArray(const MmapArray&) = delete;
   MmapArray& operator=(const MmapArray&) = delete;
      private:
        // Start of the mapped region
        char *mmap_addr;
        // Total length of the array
        size_t mmap_length;

        // Position where the first byte of the file is mapped
        char *file_mmap_addr;
        // Length of the array after the start of the file map
        size_t file_mmap_length;

      public:
        char *begin();

        MmapArray() = default;
        ~MmapArray();
        int initMmap(const std::string &filename, size_t length, size_t suffix_padding);
};

#endif  // MMAPARRAY_H_
