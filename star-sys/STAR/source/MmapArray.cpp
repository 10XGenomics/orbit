# include "MmapArray.h"

#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

MmapArray::MmapArray() {
    prefix_mmap=NULL;
    file_mmap=NULL;
    file_mmap_length=0;
};

#define PAGE (1<<12)

uint round_up_to_page(int v) {
    return v + (PAGE - 1) & !(PAGE - 1);
}

int MmapArray::makeMmap(string filename, uint length, uint suffix_padding) {

    int fd = open(filename.c_str(), O_RDONLY);

    file_mmap_length = round_up_to_page(length + suffix_padding);

    // allocate the size we need, plus a preceding page
    void* addr1 = mmap(NULL, file_mmap_length + PAGE, PROT_READ, MAP_ANONYMOUS|MAP_PRIVATE|MAP_NORESERVE, -1, 0);
    if (addr1 == (void*) -1) {
        return -1;
    }

    // map the file starting after the first page
    void* addr2 = mmap((void*) ((char*)addr1 + PAGE), file_mmap_length, PROT_READ, MAP_FIXED|MAP_PRIVATE, fd, 0);
    if (addr2 == (void*) -1) {
        return -1;
    }

    assert((char*)addr2 == (char*)addr1 + PAGE);

    // madivse to let the kernel know we want the whole file in memory
    int res = posix_madvise((char*)addr2 + PAGE, file_mmap_length, POSIX_MADV_WILLNEED);
    if (res) {
        return -1;
    }

    // save pointers
    file_mmap = addr1;
    prefix_mmap = addr2;
    
    return 0;
};

MmapArray::~MmapArray() {
    if (fd) {
        close(fd);
    }

    if (file_mmap) {
        munmap(file_mmap, (size_t) file_mmap_length);
        file_mmap=NULL;
    }

    if (prefix_mmap) {
        munmap(prefix_mmap, PAGE);
        prefix_mmap=NULL;
    }
};
