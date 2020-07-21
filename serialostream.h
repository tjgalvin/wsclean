#ifndef SERIAL_OSTREAM_H
#define SERIAL_OSTREAM_H

#include <aocommon/uvector.h>

#include <stdint.h>

#include <complex>
#include <memory>
#include <vector>

class SerialOStream {
 public:
  typedef std::size_t size_t;

  SerialOStream() {}

  size_t size() const { return _buffer.size(); }
  const unsigned char* data() const { return _buffer.data(); }
  std::string ToString() const {
    return std::string(reinterpret_cast<const char*>(_buffer.data()), size());
  }

  unsigned char* Chunk(size_t size) {
    _buffer.resize(_buffer.size() + size);
    return &*(_buffer.end() - size);
  }

  template <typename T>
  SerialOStream& UInt64(T value) {
    return write<uint64_t>(value);
  }

  template <typename T>
  SerialOStream& UInt32(T value) {
    return write<uint32_t>(value);
  }

  template <typename T>
  SerialOStream& UInt16(T value) {
    return write<uint16_t>(value);
  }

  template <typename T>
  SerialOStream& UInt8(T value) {
    return write<uint8_t>(value);
  }

  SerialOStream& Bool(bool value) { return UInt8(value ? 1 : 0); }

  SerialOStream& Float(float value) { return write(value); }

  SerialOStream& Double(double value) { return write(value); }

  SerialOStream& LDouble(long double value) { return write(value); }

  SerialOStream& CFloat(std::complex<float> value) { return write(value); }

  SerialOStream& CDouble(std::complex<double> value) { return write(value); }

  SerialOStream& CLDouble(std::complex<long double> value) {
    return write(value);
  }

  SerialOStream& String(const std::string& str) {
    UInt64(str.size());
    unsigned char* block = Chunk(str.size());
    std::copy(str.begin(), str.end(), block);
    return *this;
  }

  SerialOStream& VectorUInt64(const std::vector<uint64_t>& values) {
    // Specialization of VectorUInt64: if T==uint64_t, we don't
    // have to do conversion, so the memory block can be copied directly
    return writeVector(values);
  }

  template <typename T>
  SerialOStream& VectorUInt64(const std::vector<T>& values) {
    UInt64(values.size());
    for (const T val : values) UInt64(val);
    return *this;
  }

  SerialOStream& VectorFloat(const std::vector<float>& values) {
    return writeVector(values);
  }

  SerialOStream& VectorDouble(const std::vector<double>& values) {
    return writeVector(values);
  }

  SerialOStream& VectorCFloat(const std::vector<std::complex<float>>& values) {
    return writeVector(values);
  }

  SerialOStream& VectorCDouble(
      const std::vector<std::complex<double>>& values) {
    return writeVector(values);
  }

  template <typename T>
  SerialOStream& Ptr(const std::unique_ptr<T>& ptr) {
    return writePtr(ptr);
  }

  template <typename T>
  SerialOStream& Ptr(const std::shared_ptr<T>& ptr) {
    return writePtr(ptr);
  }

  template <typename T>
  SerialOStream& Ptr(const T* ptr) {
    return writePtr(ptr);
  }

 private:
  friend class SerialIStream;

  template <typename T>
  SerialOStream& write(T value) {
    *reinterpret_cast<T*>(Chunk(sizeof(T))) = value;
    return *this;
  }

  template <typename T>
  SerialOStream& writeVector(const std::vector<T>& values) {
    uint64_t size = values.size();
    UInt64(size);
    size_t n = values.size() * sizeof(T);
    unsigned char* block = Chunk(n);
    const unsigned char* valuePtr =
        reinterpret_cast<const unsigned char*>(values.data());
    std::copy_n(valuePtr, n, block);
    return *this;
  }

  template <typename PtrT>
  SerialOStream& writePtr(const PtrT& ptr) {
    if (ptr) {
      Bool(true);
      ptr->Serialize(*this);
    } else {
      Bool(false);
    }
    return *this;
  }

  aocommon::UVector<unsigned char> _buffer;
};

#endif
