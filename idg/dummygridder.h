#ifndef DUMMY_GRIDDER_H
#define DUMMY_GRIDDER_H

#include <complex>
#include <vector>

#include "interface.h"

namespace wsclean {

class DummyGridder : public HighLevelGridderInterface {
 public:
  virtual ~DummyGridder() {}

  virtual void set_frequencies(const double* frequency_list,
                               size_t channel_count) {
    channel_count_ = channel_count_;
  }

  virtual void set_stations(const size_t n_stations) {}

  virtual void set_kernel(size_t kernel_size, const double* kernel) {}

  virtual void start_w_layer(double layer_w_in_lambda) {}

  virtual void finish_w_layer() {}

  virtual void start_aterm(const std::complex<double>* aterm) {}

  virtual void finish_aterm() {}

  virtual void set_grid(std::complex<double>* grid) {}

  virtual void grid_visibility(
      const std::complex<float>* visibility,  // size CH x PL
      const double* uvw_in_m, size_t antenna1, size_t antenna2,
      size_t time_index) {}

  virtual void transform_grid_after_gridding() {}

  virtual void transform_grid_before_sampling() {}

  virtual void queue_visibility_sampling(const double* uvw_in_m,
                                         size_t antenna1, size_t antenna2,
                                         size_t time_index, size_t row_id,
                                         bool& is_buffer_full) {
    row_ids_.push_back(row_id);
  }

  virtual void finish_sampled_visibilities() {}

  virtual void get_sampled_visibilities(size_t index, std::complex<float>* data,
                                        size_t& rowID) const {
    for (size_t i = 0; i != channel_count_ * 4; ++i) data[i] = 1.0;
    rowID = row_ids_[index];
  }

  virtual size_t get_sampling_buffer_size() const { return row_ids_.size(); }

 private:
  std::vector<size_t> row_ids_;
  size_t channel_count_;
};

}  // namespace wsclean

#endif
