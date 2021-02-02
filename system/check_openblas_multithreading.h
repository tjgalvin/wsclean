#include <cstddef>
#include <dlfcn.h>
#include <stdexcept>

// Detect linkage to a multi-threaded version of OpenBLAS
// OpenBLAS multi-threading interfers with multi-threading in wsclean and OpenMP
void check_openblas_multithreading() {
  // Ask the dynamic linker to lookup the openblas_set_num_threads function
  void (*openblas_set_num_threads)(int) = reinterpret_cast<void (*)(int)>(
      dlsym(RTLD_DEFAULT, "openblas_set_num_threads"));

  // If openblas_set_num_threads is present, the executable is linked to a
  // multithreaded version of OpenBLAS
  if (openblas_set_num_threads) {
    // Read the OPENBLAS_NUM_THREADS environment variable
    int openblas_num_threads = 0;
    char *openblas_num_threads_env_var = getenv("OPENBLAS_NUM_THREADS");
    if (openblas_num_threads_env_var != nullptr) {
      openblas_num_threads = atoi(openblas_num_threads_env_var);
    }

    if (openblas_num_threads != 1) {
      // TODO: Fix the problem by calling openblas_set_num_threads(1), resetting
      // thread affinity. aocommon::NCPU will then return the correct value
      // again instead of 1. Then detect whether OpenMP is used (IDG uses
      // OpenMP) and then somehow reinitialize OpenMP such that not all its
      // threads are bound to CPU 0. But for now throw an error and ask the user
      // to fix the problem

      throw std::runtime_error(
          "WSClean has been linked to a multi-threaded version of OpenBLAS. \n"
          "OpenBLAS multi-threading interfers with other multi-threaded parts "
          "of WCSlean\n"
          "This has a severe impact on performance.\n"
          "Please disable OpenBLAS multi-threading by setting the environment "
          "variable OPENBLAS_NUM_THREADS to 1.\n"
          "Use either\n"
          "  setenv OPENBLAS_NUM_THREADS 1\n"
          "for csh like shells, or\n"
          "  export OPENBLAS_NUM_THREADS=1\n"
          "for bash like shells\n");

      // Alternatively the main executable wsclean could be renamed to
      // wsclean-bin And wsclean could be a wrapper shell script, that sets
      // OPENBLAS_NUM_THREAD to 1 before calling wsclean-bin
    }
  }
}
