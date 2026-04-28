#ifndef WSCLEAN_MSPROVIDERS_REORDERING_ROW_WRITER_H_
#define WSCLEAN_MSPROVIDERS_REORDERING_ROW_WRITER_H_

#include <casacore/casa/Arrays/Array.h>

namespace wsclean {

class RowWriter {
 public:
  virtual ~RowWriter() = default;

  virtual void Write(size_t row_nr, const casacore::Array<std::complex<float>>& data) = 0;
};

}  // namespace wsclean

#endif
