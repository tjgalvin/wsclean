#ifndef WSCLEAN_DECONVOLUTION_SUBIMAGELOGSET_H_
#define WSCLEAN_DECONVOLUTION_SUBIMAGELOGSET_H_

#include "controllablelog.h"

/**
 * @brief Thread safe logger for the set of subimages in the deconvolution.
 */
class SubImageLogSet {
 public:
  SubImageLogSet() : n_horizontal_(0), n_vertical_(0) {}

  void Initialize(size_t n_horizontal, size_t n_vertical) {
    std::lock_guard<std::mutex> lock(output_mutex_);
    n_horizontal_ = n_horizontal;
    n_vertical_ = n_vertical;
    const size_t n = n_horizontal * n_vertical;
    logs_.clear();
    logs_.reserve(n);
    for (size_t i = 0; i != n; ++i) {
      logs_.emplace_back(&output_mutex_);
      logs_[i].SetTag("P" + std::to_string(i) + ": ");
      logs_[i].Mute(true);
      logs_[i].Activate(false);
    }
  }

  ControllableLog& operator[](size_t index) { return logs_[index]; }

  void Activate(size_t index) {
    std::lock_guard<std::mutex> o_lock(output_mutex_);
    if (!logs_[index].IsActive()) {
      logs_[index].Activate(true);

      UnmuteMostCentral();
    }
  }

  void Deactivate(size_t index) {
    std::lock_guard<std::mutex> o_lock(output_mutex_);
    if (logs_[index].IsActive()) {
      logs_[index].Mute(true);
      logs_[index].SetOutputOnce(std::string());
      logs_[index].Activate(false);

      UnmuteMostCentral();
    }
  }

 private:
  void UnmuteMostCentral() {
    size_t unmuted_log = logs_.size();
    for (size_t i = 0; i != logs_.size(); ++i) {
      if (!logs_[i].IsMuted()) unmuted_log = i;
      logs_[i].SetTag("P" + std::to_string(i) + ": ");
      logs_[i].Mute(true);
    }

    // Find an active subimage that is as close as possible to
    // the centre (since these often contain the most interesting info)
    bool found = false;
    size_t si_x = 0;
    size_t si_y = 0;
    size_t si_d = 0;
    for (size_t y = 0; y != n_vertical_; ++y) {
      for (size_t x = 0; x != n_horizontal_; ++x) {
        const size_t index = y * n_horizontal_ + x;
        if (logs_[index].IsActive()) {
          const int dx = int(x) - int(n_horizontal_) / 2;
          const int dy = int(y) - int(n_vertical_) / 2;
          const size_t dist_squared = dx * dx + dy * dy;
          if (!found || dist_squared < si_d) {
            si_x = x;
            si_y = y;
            si_d = dist_squared;
            found = true;
          }
        }
      }
    }

    if (found) {
      const size_t index = si_y * n_horizontal_ + si_x;
      logs_[index].Mute(false);
      logs_[index].SetTag(" ");

      if (index != unmuted_log) {
        logs_[index].SetOutputOnce("Switching to output of subimage " +
                                   std::to_string(index) + "\n");
      }
    }
  }

  std::mutex output_mutex_;
  std::vector<ControllableLog> logs_;
  size_t n_horizontal_;
  size_t n_vertical_;
};

#endif