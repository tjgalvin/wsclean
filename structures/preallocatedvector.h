#ifndef WSCLEAN_STRUCTURES_PREALLOCATED_VECTOR_H_
#define WSCLEAN_STRUCTURES_PREALLOCATED_VECTOR_H_

#include <cassert>
#include <memory>
#include <utility>

/**
 * @class PreallocatedVector
 * Minimal vector-like container with preallocated fixed capacity that prevents reallocation
 *
 * This container allocates its full capacity upfront via reserve() and never reallocates,
 * making the underlying buffer address stable for concurrent access patterns when the caller
 * uses external synchronization to ensure proper memory visibility and ordering.
 *
 * Properties of this class:
 * Element access and size are not synchronized by the class.
 * Element access and a single PushBack() operation can be done in parallel as long as they do not
 * access the same (existing / newly added) elements. No reallocation after reserve(): the buffer
 * address remains constant after the initial reservation. The vector can't grow past the allocated
 * capacity (enforced in debug mode).
 *
 * This class is similar to std::vector, but std::vector does not guarantee that pushing back
 * elements at the same time as reading other values from the beginning of the vector is thread safe
 * (and in fact, in debug mode the bounds checking causes a race condition). Similarly,
 * boost::static_vector is also not thread safe in the same scenario.
 *
 * Example Usage:
 *
 * ```
 * PreallocatedVector<BeamResponse> responses;
 * responses.reserve(1000);  // Preallocate capacity - must be done before any insertions
 *
 * // Phase 1: Single thread appends data
 * for (...) { responses.PushBack(data); }
 *
 * {
 *   // Phase 2: Multiple threads modify sections (with mutex protection)
 *   {
 *       std::lock_guard<std::mutex> lock(section_mutex);
 *       for (...)
 *       {
 *           responses[index] = modified_data;  // Modify within protected section
 *       }
 *   }
 *   // Phase 3: Multiple threads read in parallel (after modifications complete)
 *   for (...)
 *   {
 *       const auto& item = responses[index];  // Safe read-only access
 *   }
 * }
 * ```
 */
template <typename T>
class PreallocatedVector {
 public:
  constexpr PreallocatedVector() noexcept = default;
  PreallocatedVector(const PreallocatedVector&) = delete;
  PreallocatedVector& operator=(const PreallocatedVector&) = delete;
  PreallocatedVector(PreallocatedVector&&) = delete;
  PreallocatedVector& operator=(PreallocatedVector&&) = delete;

  /**
   * Reserves capacity for the buffer
   *
   * Requirements:
   * - Must be called before any insertions
   * - Must not be called more than once
   */
  void Reserve(size_t new_capacity) {
    assert(size_ == 0 && "Reserve() must be called before any insertions");
    assert(capacity_ == 0 && "Reserve() may not be called once capacity is already set");
    data_ = std::make_unique<T[]>(new_capacity);
    capacity_ = new_capacity;
  }

  /**
   * Adds an element
   *
   * Requirements:
   * - Caller must ensure proper synchronization with any concurrent access
   * - Must not exceed reserved capacity
   */
  void PushBack(const T& value) {
    assert(size_ < capacity_ &&
           "PreallocatedVector capacity exceeded - increase Reserve() allocation");
    data_[size_++] = value;
  }

  /**
   * Access element for modification by index
   *
   * Requirements:
   * - Index must be within bounds (i < size())
   * - Caller must ensure concurrent modifications to the same index are impossible
   *   (e.g. by using mutexes in the calling code to protect sections)
   */
  T& operator[](size_t i) {
    assert(i < capacity_ && "Index out of bounds");
    return data_[i];
  }

  /**
   * Constant version of `operator[](size_t i)` (see above).
   */
  const T& operator[](size_t i) const {
    assert(i < capacity_ && "Index out of bounds");
    return data_[i];
  }

  /**
   * Returns current size
   *
   * Note: In concurrent scenarios, caller should synchronize this call.
   */
  size_t Size() const { return size_; }

  size_t Capacity() const { return capacity_; }

  bool Empty() const { return size_ == 0; }

 private:
  std::unique_ptr<T[]> data_ = nullptr;
  size_t size_ = 0;
  size_t capacity_ = 0;
};

#endif  // WSCLEAN_STRUCTURES_PREALLOCATED_VECTOR_H_
