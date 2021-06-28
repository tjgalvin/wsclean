#ifndef NOISE_MS_ROW_PROVIDER_H
#define NOISE_MS_ROW_PROVIDER_H

#include "directmsrowprovider.h"

#include <map>
#include <random>
#include <string>

class NoiseMSRowProvider : public DirectMSRowProvider {
 public:
  NoiseMSRowProvider(const string& msPath, const MSSelection& selection,
                     const std::map<size_t, size_t>& selectedDataDescIds,
                     const std::string& dataColumnName, bool requireModel);

  void SetNoiseLevel(double noiseStdDevJy);

  void SetNoiseBaselineFile(const std::string& filename);

  virtual void ReadData(DataArray& data, FlagArray& flags, WeightArray& weights,
                        double& u, double& v, double& w, uint32_t& dataDescId,
                        uint32_t& antenna1, uint32_t& antenna2,
                        uint32_t& fieldId, double& time) final override;

  class NoiseMap {
   public:
    NoiseMap() = default;
    NoiseMap(std::istream& stream);
    bool Empty() const { return _map.empty(); }
    float GetNoiseValue(size_t antenna1, size_t antenna2) const;

   private:
    // Maps an antenna1, antenna2 pair to the stddev level
    std::map<std::pair<size_t, size_t>, float> _map;
  };

 private:
  std::mt19937 _rng;
  std::normal_distribution<float> _distribution;
  NoiseMap _noiseMap;
};

#endif
