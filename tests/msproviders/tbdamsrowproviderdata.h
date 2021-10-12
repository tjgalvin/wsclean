#ifndef msproviders_tbdamsrowproviderdata_h
#define msproviders_tbdamsrowproviderdata_h

#include <array>
#include <cinttypes>

namespace MwaBdaMockMs {
enum Columns : std::size_t { kTime = 0, kAntenna1 = 1, kAntenna2 = 2 };
extern const std::array<std::array<uint64_t, 3>, 21> kMs;
}  // namespace MwaBdaMockMs

#endif  // msproviders_tbdamsrowproviderdata_h
