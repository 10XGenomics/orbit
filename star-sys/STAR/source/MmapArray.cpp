# include "MmapArray.h"

#include <string>
#include <iostream>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#define PAGE ((size_t) (1<<12))

size_t round_up_to_page(size_t v) {
    return (v + (PAGE - 1)) & ~(PAGE-1);
}

// mmap `filename` as a mostly read-only map, allowing for limited writing before an after the file.
// Maps at least `length` bytes of the file, although additional bytes from the file may be mapped.
// Ensures that there is one writeable page before address where the file is mapped.
// Ensures that there is at least `suffix_padding` writable bytes beyond the requested mapping length.
// madvise()s the file mapping to instruct the kernel to begin reading the data.
int MmapArray::initMmap(const std::string &filename, size_t length, size_t suffix_padding) {    

    int fd = open(filename.c_str(), O_RDONLY|O_CLOEXEC);
    if (fd == -1) {
        return -1;
    }

    file_mmap_length = round_up_to_page(length + suffix_padding);
    mmap_length = file_mmap_length + PAGE;

    // allocate the size we need, plus a preceding page
    void* addr1 = mmap(NULL, mmap_length, PROT_READ|PROT_WRITE, MAP_ANONYMOUS|MAP_PRIVATE|MAP_NORESERVE, -1, 0);

    if (addr1 == (void*) -1) {
        int err = errno;
        std::cerr << "reference mmap failed: " << strerror(err) << " on fd: " << fd;
        std::cerr << " len pre: " << length << " len_round: " << file_mmap_length;
        std::cerr << " addr1: " << addr1 << "\n";
        close(fd);
        return -1;
    }

    // map the file starting after the first page
    void* addr2 = mmap((void*) ((char*)addr1 + PAGE), file_mmap_length, PROT_READ, MAP_FIXED|MAP_PRIVATE, fd, 0);
    
    // make sure the mapping went to the address we requested.
    if (addr2 == (void*) -1 || (char*)addr2 != (char*)addr1 + PAGE) {
        int err = errno;
        std::cerr << "reference mmap failed: " << strerror(err) << " on fd: " << fd;
        std::cerr << " : addr2: " << addr2 << "\n";
        munmap(mmap_addr, mmap_length);
        close(fd);
        return -1;
    }

    // make the final page writeable. On Mac at least, it's much faster if the 'main' mmap above is read only.
    void* addr3 = mmap((void*) ((char*)addr1 + file_mmap_length), PAGE, PROT_READ|PROT_WRITE, MAP_FIXED|MAP_PRIVATE, fd, file_mmap_length - PAGE);

    if (addr3 == (void*) -1 || (char*)addr3 != (char*)addr1 + file_mmap_length) {
        int err = errno;
        std::cerr << "reference mmap failed: " << strerror(err) << " on fd: " << fd;
        std::cerr << "reference mmap failed: addr3: " << addr3 << "\n";
        munmap(mmap_addr, mmap_length);
        close(fd);
        return -1;
    }

    // don't need the file handle anymore
    close(fd);

    // madivse to let the kernel know we want the whole file in memory
    int res = posix_madvise((char*)addr2, file_mmap_length, POSIX_MADV_WILLNEED);
    
    if (res) {
        std::cerr << "madivse on reference mmap failred: res: " << res << "\n";
    }

    // save pointers
    mmap_addr = (char*) addr1;
    file_mmap_addr = (char*) addr2;
    
    return 0;
};

char* MmapArray::begin() {
    return file_mmap_addr;
}

MmapArray::~MmapArray() {
    // this unmaps everything in one shot
    if (mmap_addr) {
        munmap(mmap_addr, mmap_length);
        mmap_addr=nullptr;
    }
};
