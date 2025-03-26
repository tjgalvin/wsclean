#include "wgridder_implementation.h"

template class wsclean::WGridder<float>;

template <typename NumT>
template <typename Tms>
void wsclean::WGridder<NumT>::AddInversionMs(size_t n_rows, const double *uvw,
                                             const ducc0::cmav<double, 1> &freq,
                                             Tms &ms) {
  AddInversionMsImplementation(n_rows, uvw, freq, ms);
}
template void wsclean::WGridder<float>::AddInversionMs<
    wsclean::VisibilityCallbackBuffer<std::complex<float>> const>(
    size_t n_rows, const double *uvw, const ducc0::cmav<double, 1> &freq,
    wsclean::VisibilityCallbackBuffer<std::complex<float>> const &ms);
