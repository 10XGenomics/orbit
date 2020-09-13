# include "MmapArray.h"

#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

MmapArray::MmapArray() {
    prefix_mmap=NULL;
    file_mmap=NULL;
    file_mmap_length=0;
    fd=0;
};

#define PAGE ((uint) (1<<12))

uint round_up_to_page(uint v) {
    return v + (PAGE - 1) & ~(PAGE-1);
}

int MmapArray::makeMmap(string filename, uint length, uint suffix_padding) {

    fd = open(filename.c_str(), O_RDONLY);
    if (fd == -1) {
        return -1;
    }

    file_mmap_length = round_up_to_page(length + suffix_padding);
    
    // allocate the size we need, plus a preceding page
    void* addr1 = mmap(NULL, file_mmap_length + PAGE, PROT_READ|PROT_WRITE, MAP_ANONYMOUS|MAP_PRIVATE|MAP_NORESERVE, -1, 0);

    if (addr1 == (void*) -1) {
        std::cerr << "fd: " << fd;
        std::cerr << "len pre: " << length << "len_round: " << file_mmap_length << "\n";
        std::cerr << "addr1: " << addr1;
        return -1;
    }

    // map the file starting after the first page
    void* addr2 = mmap((void*) ((char*)addr1 + PAGE), file_mmap_length, PROT_READ, MAP_FIXED|MAP_PRIVATE, fd, 0);
    
    if (addr2 == (void*) -1) {
        std::cerr << "addr2: " << addr2;
        return -1;
    }

    assert((char*)addr2 == (char*)addr1 + PAGE);

    // make the final page writeable. On Mac at least, it's much faster if the 'main' mmap above is read only.
    void* addr3 = mmap((void*) ((char*)addr1 + file_mmap_length), PAGE, PROT_READ|PROT_WRITE, MAP_FIXED|MAP_PRIVATE, fd, 0);
    
    if (addr3 == (void*) -1) {
        std::cerr << "addr3: " << addr3;
        return -1;
    }

    assert((char*)addr3 == (char*)addr1 + file_mmap_length);


    // madivse to let the kernel know we want the whole file in memory
    int res = posix_madvise((char*)addr2, file_mmap_length, POSIX_MADV_WILLNEED);
    
    if (res) {
        std::cerr << "res: " << res;
        return -1;
    }

    // save pointers
    file_mmap = (char*) addr2;
    prefix_mmap = (char*) addr1;
    
    return 0;
};

MmapArray::~MmapArray() {

    if (file_mmap) {
        munmap(file_mmap, (size_t) file_mmap_length);
        file_mmap=NULL;
    }

    if (prefix_mmap) {
        munmap(prefix_mmap, PAGE);
        prefix_mmap=NULL;
    }

    if (last_page_mmap) {
        munmap(last_page_mmap, PAGE);
        last_page_mmap=NULL;
    }

    if (fd) {
        close(fd);
    }
};
