#include <boost/test/unit_test.hpp>

namespace test {

/**
 * Helper class for using using unique pointers in tests.
 *
 * Use case: When using unique pointers, tests often also store raw pointers
 * for accessing objects after std::moving their unique pointer.
 *
 * This class supports this common pattern: Instead of storing both a unique
 * pointer and a raw pointer, tests can use a single test::UniquePtr variable.
 *
 * The perfect forwarding constructor also simplifies initialization:
 * "std::unique_ptr<X> x(new X(...));" becomes: "test::UniquePtr<X> x(...);"
 *
 * @tparam T Object type for the unique pointer.
 */
template <class T>
class UniquePtr {
 public:
  /**
   * Constructor. It uses perfect forwarding for calling T's constructor.
   */
  template <typename... Args>
  UniquePtr(Args&&... args)
      : _unique(new T(std::forward<Args>(args)...)), _raw(_unique.get()) {}

  /**
   * Extract the unique_ptr from the helper class.
   * The test code may only call this function once.
   */
  std::unique_ptr<T> take() {
    BOOST_TEST_REQUIRE(static_cast<bool>(_unique));
    return std::move(_unique);
  }

  T& operator*() { return *_raw; }
  T* operator->() { return _raw; }
  T* get() { return _raw; }

 private:
  std::unique_ptr<T> _unique;
  T* const _raw;
};

}  // namespace test
