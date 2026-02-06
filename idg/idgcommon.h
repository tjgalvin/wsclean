#ifndef IDG_MS_GRIDDER_COMMON_H
#define IDG_MS_GRIDDER_COMMON_H

#include <idg-api.h>

#include "../main/settings.h"
#include "../gridding/msgridderdata.h"

namespace wsclean {

struct IDGInversionRow final : public MsGridderData::InversionRow {
  size_t antenna1;
  size_t antenna2;
  size_t time_index;
};

struct IDGPredictionRow {
  double uvw[3];

  size_t antenna1;
  size_t antenna2;
  size_t time_index;
  size_t row_id;
};

idg::api::Type GetIdgType(const Settings& settings);

}  // namespace wsclean

#endif
