#ifndef WSCLEAN_GRIDDING_H5_SOLUTION_DATA_H_
#define WSCLEAN_GRIDDING_H5_SOLUTION_DATA_H_

#include <vector>

#include <schaapcommon/h5parm/h5parm.h>
#include <schaapcommon/h5parm/jonesparameters.h>
#include <schaapcommon/h5parm/soltab.h>

#include "../main/settings.h"

namespace wsclean {

class GriddingTaskManager;

class H5SolutionData {
 public:
  H5SolutionData(const Settings& settings);

  H5SolutionData(const H5SolutionData&) = delete;
  H5SolutionData& operator=(const H5SolutionData&) = delete;

  bool HasData() const { return !h5parms_.empty(); }

  const std::vector<schaapcommon::h5parm::H5Parm>& GetH5Parms() const { return h5parms_; }
  const std::vector<schaapcommon::h5parm::SolTab*>& GetFirstSolutions() const {
    return first_solutions_;
  }
  const std::vector<schaapcommon::h5parm::SolTab*>& GetSecondSolutions() const {
    return second_solutions_;
  }
  const std::vector<schaapcommon::h5parm::GainType>& GetGainTypes() const { return gain_types_; }

 private:
  void LoadSolutions();
  void LoadGainTypes();

  const Settings& settings_;

  /** For each solution, an H5Parm object for accessing the values. */
  std::vector<schaapcommon::h5parm::H5Parm> h5parms_;

  /**
   * For each solution, a pointer into h5parms_ for the first part.
   * If there are two solution tables, it contains the amplitudes.
   */
  std::vector<schaapcommon::h5parm::SolTab*> first_solutions_;

  /**
   * For each solution, a pointer into h5parms_ for the second (phase) part.
   * This vector remains empty if there is a single solution table.
   */
  std::vector<schaapcommon::h5parm::SolTab*> second_solutions_;

  /** For each solution, the corresponding gain type. */
  std::vector<schaapcommon::h5parm::GainType> gain_types_;
};

}  // namespace wsclean

#endif  // WSCLEAN_GRIDDING_H5_SOLUTION_DATA_H_
