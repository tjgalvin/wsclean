#include "idgcommon.h"

namespace wsclean {

idg::api::Type GetIdgType(const Settings& settings) {
  switch (settings.idgMode) {
    case Settings::IDG_GPU:
      return idg::api::Type::CUDA_GENERIC;
    case Settings::IDG_HYBRID:
      return idg::api::Type::HYBRID_CUDA_CPU_OPTIMIZED;
    case Settings::IDG_CPU:
    default:
      return idg::api::Type::CPU_OPTIMIZED;
  }
}

}  // namespace wsclean
