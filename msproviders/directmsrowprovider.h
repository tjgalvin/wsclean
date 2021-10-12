#ifndef DIRECT_MS_ROW_PROVIDER_H
#define DIRECT_MS_ROW_PROVIDER_H

#include "msrowprovider.h"

class DirectMSRowProvider : public MSRowProvider {
 public:
  DirectMSRowProvider(const string& msPath, const MSSelection& selection,
                      const std::map<size_t, size_t>& selectedDataDescIds,
                      const std::string& dataColumnName, bool requireModel)
      : MSRowProvider(msPath, selection, selectedDataDescIds, dataColumnName,
                      requireModel) {}

  explicit DirectMSRowProvider(
      const casacore::MeasurementSet& ms, const MSSelection& selection,
      const std::map<size_t, size_t>& selected_data_description_ids,
      const std::string& data_column_name, bool require_model)
      : MSRowProvider(ms, selection, selected_data_description_ids,
                      data_column_name, require_model) {}

  virtual void ReadData(DataArray& data, FlagArray& flags, WeightArray& weights,
                        double& u, double& v, double& w, uint32_t& dataDescId,
                        uint32_t& antenna1, uint32_t& antenna2,
                        uint32_t& fieldId, double& time) override;

  virtual void ReadModel(DataArray& model) final override;

 private:
};

#endif
