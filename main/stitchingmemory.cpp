#include "stitchingmemory.h"

#include <algorithm>

namespace wsclean::detail {

namespace {

inline size_t CeilDiv(size_t numerator, size_t denominator) {
  return denominator == 0 ? 0 : (numerator + denominator - 1) / denominator;
}

}  // namespace

size_t EstimateStitchJobBytes(const StitchMemoryModel& model) {
  if (model.image_width == 0 || model.image_height == 0) return 0;

  const size_t image_slot_bytes = model.image_width * model.image_height * sizeof(float);
  size_t image_count = 2;

  if (model.has_feathering && model.type != StitchImageType::kFacetPbModel) {
    image_count += 2;
  }
  if (model.apply_matrix) {
    image_count += std::max(model.polarization_count, size_t{4});
  }
  if (model.type == StitchImageType::kPSF && model.save_uv) {
    image_count += 2;
  }

  return (image_count * image_slot_bytes * 5 + 3) / 4;
}

size_t CalculateRequiredStitchMemorySlots(size_t estimated_job_bytes, size_t image_slot_bytes) {
  if (estimated_job_bytes == 0 || image_slot_bytes == 0) return 1;
  return std::max<size_t>(1, CeilDiv(estimated_job_bytes, image_slot_bytes));
}

void MemoryReservation::Release() {
  if (gate_ != nullptr) {
    gate_->Release(slots_);
    gate_ = nullptr;
    slots_ = 0;
  }
}

MemoryGate::MemoryGate(size_t total_slots)
    : total_slots_(std::max<size_t>(1, total_slots)), available_slots_(total_slots_) {}

MemoryReservation MemoryGate::Acquire(size_t slots) {
  slots = std::clamp<size_t>(slots, 1, total_slots_);
  std::unique_lock<std::mutex> lock(mutex_);
  condition_.wait(lock, [&]() { return available_slots_ >= slots; });
  available_slots_ -= slots;
  return MemoryReservation(this, slots);
}

size_t MemoryGate::AvailableSlots() const {
  std::lock_guard<std::mutex> lock(mutex_);
  return available_slots_;
}

void MemoryGate::Release(size_t slots) {
  {
    std::lock_guard<std::mutex> lock(mutex_);
    available_slots_ = std::min(total_slots_, available_slots_ + slots);
  }
  condition_.notify_all();
}

}  // namespace wsclean::detail
