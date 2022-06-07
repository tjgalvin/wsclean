#ifndef MAPPED_FILE_H_
#define MAPPED_FILE_H_

#include <aocommon/system.h>

#include <stdexcept>
#include <string>

#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <errno.h>
#include <fcntl.h>

/**
 * A memory-mapped buffer that is mapped to a file on disk using mmap().
 */
class MappedFile {
 public:
  MappedFile() : reserved_size_(0) {}
  MappedFile(const MappedFile&) = delete;
  MappedFile(MappedFile&& source)
      : reserved_size_(source.reserved_size_),
        memory_map_(source.memory_map_),
        file_descriptor_(source.file_descriptor_) {
    source.reserved_size_ = 0;
    source.memory_map_ = nullptr;
    source.file_descriptor_ = -1;
  }

  MappedFile& operator=(MappedFile&& rhs) {
    Swap(*this, rhs);
    return *this;
  }
  MappedFile& operator=(const MappedFile& rhs) = delete;

  /**
   * Memory map a file. The file should exist prior to this call,
   * and should contain data.
   */
  MappedFile(const std::string& filename, size_t reserved_size)
      : reserved_size_(reserved_size) {
    file_descriptor_ = open(filename.c_str(), O_RDWR);
    if (file_descriptor_ == -1)
      throw std::runtime_error("Error opening file '" + filename + "'");
    if (reserved_size_ != 0) {
#ifdef MAP_NORESERVE
      memory_map_ = reinterpret_cast<char*>(
          mmap(nullptr, reserved_size_, PROT_READ | PROT_WRITE,
               MAP_SHARED | MAP_NORESERVE, file_descriptor_, 0));
#else
      memory_map_ = reinterpret_cast<char*>(
          mmap(nullptr, reserved_size_, PROT_READ | PROT_WRITE, MAP_SHARED,
               file_descriptor_, 0));
#endif
      if (memory_map_ == MAP_FAILED) {
        memory_map_ = nullptr;
        close(file_descriptor_);
        const std::string msg = aocommon::system::GetErrorString(errno);
        throw std::runtime_error(
            std::string("Error creating memory map for file '" + filename +
                        "': "
                        "mmap() returned MAP_FAILED with error message: ") +
            msg);
      }
    }
  }

  ~MappedFile() {
    if (memory_map_ != nullptr) {
      if (reserved_size_ != 0) munmap(memory_map_, reserved_size_);
    }
    if (file_descriptor_ != -1) {
      close(file_descriptor_);
    }
  }

  char* Data() { return memory_map_; }
  const char* Data() const { return memory_map_; }

 private:
  friend void Swap(MappedFile& lhs, MappedFile& rhs) {
    std::swap(lhs.reserved_size_, rhs.reserved_size_);
    std::swap(lhs.memory_map_, rhs.memory_map_);
    std::swap(lhs.file_descriptor_, rhs.file_descriptor_);
  }

  size_t reserved_size_;
  char* memory_map_ = nullptr;
  int file_descriptor_ = -1;
};

#endif
