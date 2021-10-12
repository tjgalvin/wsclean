#include "tbdamsrowproviderdata.h"

namespace MwaBdaMockMs {
// clang-format off
// The data is generated using the following TaQL query:
// 'select int(TIME),ANTENNA1,ANTENNA2 from MWA_BDA_MOCK.ms where ANTENNA1 != ANTENNA2'
// clang-format on
// The MWA_BDA_MOCK.ms data set is downloaded during the build process of
// WSClean. Since the processing of the MS in WSClean always skips the
// autocorrelation these entires are filtered out by the query.
const std::array<std::array<uint64_t, 3>, 21> kMs{
    {{4875418082, 97, 126},  {4875418084, 1, 97},    {4875418084, 2, 97},
     {4875418084, 87, 126},  {4875418086, 1, 78},    {4875418086, 1, 87},
     {4875418086, 1, 122},   {4875418086, 2, 78},    {4875418086, 2, 87},
     {4875418086, 2, 122},   {4875418086, 64, 126},  {4875418086, 81, 126},
     {4875418086, 83, 126},  {4875418086, 84, 126},  {4875418090, 97, 126},
     {4875418086, 109, 126}, {4875418086, 110, 126}, {4875418086, 111, 126},
     {4875418092, 1, 97},    {4875418092, 2, 97},    {4875418092, 87, 126}}};

}  // namespace MwaBdaMockMs
