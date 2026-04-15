#ifndef STITCHING_MEMORY_H
#define STITCHING_MEMORY_H

#include <condition_variable>
#include <cstddef>
#include <mutex>

namespace wsclean {

enum class StitchImageType { kPSF, kDirty, kFacetPbModel, kOther };

namespace detail {

/**
 * Small model used to estimate the temporary memory footprint of a stitch job.
 */
struct StitchMemoryModel {
  size_t image_width = 0;
  size_t image_height = 0;
  size_t polarization_count = 0;
  StitchImageType type = StitchImageType::kOther;
  bool apply_matrix = false;
  bool has_feathering = false;
  bool save_uv = false;
};

size_t EstimateStitchJobBytes(const StitchMemoryModel& model);
size_t CalculateRequiredStitchMemorySlots(size_t estimated_job_bytes, size_t image_slot_bytes);

class MemoryGate;

/**
 * RAII handle that returns reserved slots to a MemoryGate on destruction.
 */
class MemoryReservation {
 public:
  MemoryReservation() = default;
  MemoryReservation(MemoryGate* gate, size_t slots) : gate_(gate), slots_(slots) {}
  MemoryReservation(const MemoryReservation&) = delete;
  MemoryReservation& operator=(const MemoryReservation&) = delete;
  MemoryReservation(MemoryReservation&& other) noexcept : gate_(other.gate_), slots_(other.slots_) {
    other.gate_ = nullptr;
    other.slots_ = 0;
  }
  MemoryReservation& operator=(MemoryReservation&& other) noexcept {
    if (this != &other) {
      Release();
      gate_ = other.gate_;
      slots_ = other.slots_;
      other.gate_ = nullptr;
      other.slots_ = 0;
    }
    return *this;
  }
  ~MemoryReservation() noexcept { Release(); }

  size_t Slots() const { return slots_; }
  void Release();

 private:
  MemoryGate* gate_ = nullptr;
  size_t slots_ = 0;
};

/**
 * Counting gate used to throttle concurrent memory-heavy stitch jobs.
 */
class MemoryGate {
 public:
  explicit MemoryGate(size_t total_slots);

  [[nodiscard]] MemoryReservation Acquire(size_t slots);
  size_t TotalSlots() const { return total_slots_; }
  size_t AvailableSlots() const;

 private:
  friend class MemoryReservation;

  void Release(size_t slots);

  size_t total_slots_;
  mutable std::mutex mutex_;
  std::condition_variable condition_;
  size_t available_slots_;
};

}  // namespace detail

}  // namespace wsclean

#endif
