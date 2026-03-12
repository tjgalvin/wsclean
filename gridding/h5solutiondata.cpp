#include "h5solutiondata.h"

#include <vector>

#include <schaapcommon/h5parm/h5parm.h>
#include <schaapcommon/h5parm/jonesparameters.h>
#include <schaapcommon/h5parm/soltab.h>

#include "../main/settings.h"

namespace wsclean {

H5SolutionData::H5SolutionData(const Settings& settings) : settings_(settings) {
  LoadSolutions();
  LoadGainTypes();
}

void H5SolutionData::LoadSolutions() {
  using schaapcommon::h5parm::SolTab;

  const std::vector<std::string>& solution_files = settings_.facetSolutionFiles;
  const std::vector<std::string>& solution_table_names = settings_.facetSolutionTables;

  h5parms_.reserve(solution_files.size());
  first_solutions_.reserve(solution_files.size());
  if (solution_table_names.size() == 2) second_solutions_.reserve(solution_files.size());

  for (const std::string& solution_file : solution_files) {
    schaapcommon::h5parm::H5Parm& h5parm =
        h5parms_.emplace_back(solution_file, solution_table_names);

    if (solution_table_names.size() == 1) {
      first_solutions_.emplace_back(&h5parm.GetSolTab(solution_table_names[0]));
    } else {  // Use amplitude as first solution and phase as second solution.
      assert(solution_table_names.size() == 2);
      const std::string kAmplitude = "amplitude";
      const std::string kPhase = "phase";

      const std::array<schaapcommon::h5parm::SolTab*, 2> tables{
          &h5parm.GetSolTab(solution_table_names[0]), &h5parm.GetSolTab(solution_table_names[1])};
      const std::array<std::string, 2> types{tables[0]->GetType(), tables[1]->GetType()};

      if (types[0] == kAmplitude && types[1] == kPhase) {
        first_solutions_.emplace_back(tables[0]);
        second_solutions_.emplace_back(tables[1]);
      } else if (types[0] == kPhase && types[1] == kAmplitude) {
        first_solutions_.emplace_back(tables[1]);
        second_solutions_.emplace_back(tables[0]);
      } else {
        throw std::runtime_error("WSClean expects solution tables with names '" + kAmplitude +
                                 "' and '" + kPhase + "', but received '" + types[0] + "' and '" +
                                 types[1] + "'");
      }
    }
  }
}

void H5SolutionData::LoadGainTypes() {
  using schaapcommon::h5parm::GainType;
  using schaapcommon::h5parm::SolTab;

  gain_types_.reserve(h5parms_.size());

  if (second_solutions_.empty()) {
    for (const SolTab* table : first_solutions_) {
      gain_types_.emplace_back(
          schaapcommon::h5parm::JonesParameters::H5ParmTypeStringToGainType(table->GetType()));
    }
  } else {
    assert(second_solutions_.size() == first_solutions_.size());

    const std::string kPol = "pol";

    for (size_t i = 0; i < first_solutions_.size(); ++i) {
      const size_t n_amplitude_polarizations =
          first_solutions_[i]->HasAxis(kPol) ? first_solutions_[i]->GetAxis(kPol).size : 1;
      const size_t n_phase_polarizations =
          second_solutions_[i]->HasAxis(kPol) ? second_solutions_[i]->GetAxis(kPol).size : 1;

      if (n_amplitude_polarizations == 1 && n_phase_polarizations == 1) {
        gain_types_.emplace_back(GainType::kScalarComplex);
      } else if (n_amplitude_polarizations == 2 && n_phase_polarizations == 2) {
        gain_types_.emplace_back(GainType::kDiagonalComplex);
      } else if (n_amplitude_polarizations == 4 && n_phase_polarizations == 4) {
        gain_types_.emplace_back(GainType::kFullJones);
      } else {
        throw std::runtime_error(
            "Incorrect or mismatching number of polarizations in the "
            "provided amplitude and phase soltabs. The number of polarizations "
            "should both be either 1, 2 or 4, but received " +
            std::to_string(n_amplitude_polarizations) + " for amplitude and " +
            std::to_string(n_phase_polarizations) + " for phase");
      }
    }
  }
}

}  // namespace wsclean
