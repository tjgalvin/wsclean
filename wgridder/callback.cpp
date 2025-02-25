#include "wgridder.h"

#include <complex>

#include "ducc0/wgridder/wgridder.h"

#include "../gridding/gainmode.h"

using namespace ducc0;

/*
 * This file contains the implementation of @ref VisibilityCallbackBuffer and
 * some of the methods from @ref WGridder used to
 * instantiate/call @ref VisibilityCallbackBuffer.
 * They are implemented here instead of directly in wgridder.h or wgridder.cpp
 * as they need to be explicitely instantiated multiple times with different
 * template parameters.
 * This is done by the build system compiling this file seperately for each
 * possible combination.
 * The reason for doing this is that each instantiation must compile all of
 * ducc0; which is quite a resource intensive task even when included once, but
 * in this case we have N! combinations (NumT * GainMode * NPolarizations *
 * NParms * ApplyBeam * ApplyForward * HasH5Parm).
 * Effectively hundreds of different copies of ducc0 instantiated, when they are
 * all in one file, it overwhelms the compiler, uses too much memory and
 * compiles slowly. By instantiating instead N files with a smaller subset of X
 * instantiations per file the time and memory requirements for compiling are
 * kept more manageable and the build can make better use of parallel
 * processing.
 */

namespace wsclean {

class MsGridder;

/**
 * VisibilityCallbackBuffer implements a virtual buffer replacement to the
 * `cmav` that would ordinarily be used to pass visibility data into DUCC.
 *
 * Ordinarily the `cmav` that DUCC takes would contain visibilities with facet
 * solutions pre-applied.
 * With VisibilityCallbackBuffer we instead hold in memory a buffer that does
 * not have facet solutions applied.
 * When DUCC requests from the buffer a specific visibility for a specific facet
 * the facet solution is applied "on the fly" and the required value returned.
 * Some internal caching is applied at the row level to help a bit with
 * efficiency.
 */
template <typename TVisibility, GainMode Mode, size_t NPolarizations,
          size_t NParms, bool ApplyBeam, bool ApplyForward, bool HasH5Parm,
          typename TInfo = ducc0::detail_mav::mav_info<2>>
class VisibilityCallbackBuffer : public TInfo {
 public:
  VisibilityCallbackBuffer(size_t n_rows, VisibilityCallbackData &data)
      : TInfo({n_rows, data.n_channels}),
        n_antennas_(data.n_antennas),
        n_channels_(data.n_channels),
        n_visibilities_per_row_(data.n_channels * NPolarizations),
        selected_band_(data.selected_band),
        antennas_(data.antennas),
        visibilities_(data.visibilities),
        time_offsets_(data.time_offsets),
        gridder_(data.gridder),
        parm_response_(data.parm_response) {}

  template <typename Index>
  const TVisibility raw(Index index) const {
    // Calculate offsets
    const size_t row = index / n_channels_;
    size_t channel = index % n_channels_;
    // Retrieve value for offsets
    const std::pair<size_t, size_t> &antenna_pair = antennas_[row];
    const size_t &time_offset = time_offsets_[row];
    std::complex<float> visibilities[NPolarizations];
    std::copy_n(&visibilities_[(row * n_visibilities_per_row_) +
                               (channel * NPolarizations)],
                NPolarizations, visibilities);
    // Apply correction
    gridder_->ApplySingleCorrection<Mode, NParms, ModifierBehaviour::kApply,
                                    ApplyBeam, ApplyForward, HasH5Parm>(
        parm_response_, channel, n_channels_, n_antennas_, visibilities,
        nullptr, antenna_pair.first, antenna_pair.second, time_offset, nullptr);
    internal::CollapseData<NPolarizations>(1, visibilities,
                                           gridder_->Polarization());
    return visibilities[0];
  }
  template <typename... Params>
  const TVisibility operator()(Params... params) const {
    return raw(TInfo::idx(params...));
  }

  // Turn all prefetch operations inside DUCC into null ops
  // As we return by value and are not a persistent buffer prefetching doesn't
  // make sense in this context
  template <typename Index>
  void prefetch_r(Index) const {}
  template <typename Index>
  void prefetch_w(Index) const {}
  template <typename... Params>
  void prefetch_r(Params...) const {}

 private:
  size_t n_antennas_;
  // Number of channels per row of visibilities
  size_t n_channels_;
  size_t n_visibilities_per_row_;
  const aocommon::BandData &selected_band_;
  const std::pair<size_t, size_t> *antennas_;
  const std::complex<float> *visibilities_;
  /**
   * When applying corrections sequentially a time_offset is calculated by @ref
   * CacheParmResponse() for each row, used for applying the corrections, and
   * then the time_offset for the next row calculated on top of it.
   * As the time_offset is needed when we apply the corrections, and we can't
   * compute it again here without sequentially going through every single row,
   * we have to store all of them in a buffer to be used when we apply the
   * corrections.
   */
  const size_t *time_offsets_;
  MsGridder *gridder_;
  const std::complex<float> *parm_response_;
};

template <typename NumT>
template <GainMode Mode, size_t NPolarizations, size_t NParms>
void WGridder<NumT>::CreateAndAddInversionMs4(
    bool apply_beam, bool apply_forward, bool has_h5_parm, size_t n_rows,
    const double *uvws, const ducc0::cmav<double, 1> &frequencies,
    VisibilityCallbackData &data) {
  if (apply_beam) {
    CreateAndAddInversionMs5<Mode, NPolarizations, NParms, true>(
        apply_forward, has_h5_parm, n_rows, uvws, frequencies, data);
  } else {
    CreateAndAddInversionMs5<Mode, NPolarizations, NParms, false>(
        apply_forward, has_h5_parm, n_rows, uvws, frequencies, data);
  }
}

template <typename NumT>
template <GainMode Mode, size_t NPolarizations, size_t NParms, bool ApplyBeam>
void WGridder<NumT>::CreateAndAddInversionMs5(
    bool apply_forward, bool has_h5_parm, size_t n_rows, const double *uvws,
    const ducc0::cmav<double, 1> &frequencies, VisibilityCallbackData &data) {
  if (apply_forward) {
    CreateAndAddInversionMs6<Mode, NPolarizations, NParms, ApplyBeam, true>(
        has_h5_parm, n_rows, uvws, frequencies, data);
  } else {
    CreateAndAddInversionMs6<Mode, NPolarizations, NParms, ApplyBeam, false>(
        has_h5_parm, n_rows, uvws, frequencies, data);
  }
}

template <typename NumT>
template <GainMode Mode, size_t NPolarizations, size_t NParms, bool ApplyBeam,
          bool ApplyForward>
void WGridder<NumT>::CreateAndAddInversionMs6(
    bool has_h5_parm, size_t n_rows, const double *uvws,
    const ducc0::cmav<double, 1> &frequencies, VisibilityCallbackData &data) {
  if (has_h5_parm) {
    CreateAndAddInversionMs7<Mode, NPolarizations, NParms, ApplyBeam,
                             ApplyForward, true>(n_rows, uvws, frequencies,
                                                 data);
  } else {
    CreateAndAddInversionMs7<Mode, NPolarizations, NParms, ApplyBeam,
                             ApplyForward, false>(n_rows, uvws, frequencies,
                                                  data);
  }
}

template <typename NumT>
template <GainMode Mode, size_t NPolarizations, size_t NParms, bool ApplyBeam,
          bool ApplyForward, bool HasH5Parm>
void WGridder<NumT>::CreateAndAddInversionMs7(
    size_t n_rows, const double *uvws,
    const ducc0::cmav<double, 1> &frequencies, VisibilityCallbackData &data) {
  const VisibilityCallbackBuffer<std::complex<float>, Mode, NPolarizations,
                                 NParms, ApplyBeam, ApplyForward, HasH5Parm>
      ms(n_rows, data);
  AddInversionMs(n_rows, uvws, frequencies, ms);
}

// Force explicit concrete instantiation of
// WGridder<TYPE>::CreateAndAddInversionMs2<MODE>()
// CMake will recompile this file multiple times, one for each possible
// combination. This allows us to keep each instantiation in a seperate object
// file while keeping the repeated boilerplate overhead as little as possible.
// Its desirable to have each instantiation seperate because if they are all in
// one file (as would ordinarily occur) compile times balloon out of control.

template void WGridder<DTYPE>::CreateAndAddInversionMs4<MODE, NPOL, NPARMS>(
    bool apply_beam, bool apply_forward, bool has_h5_parm, size_t n_rows,
    const double *uvws, const cmav<double, 1> &frequencies,
    VisibilityCallbackData &data);

}  // namespace wsclean
