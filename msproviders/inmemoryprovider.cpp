#include "inmemoryprovider.h"

#include "msreaders/inmemoryreader.h"

namespace wsclean {

std::unique_ptr<MSReader> InMemoryProvider::MakeReader() {
  return std::make_unique<InMemoryReader>(this, data_.rows,
                                          data_.meta_data->begin());
}

}  // namespace wsclean
