#include "noisemsrowprovider.h"

#include "../io/logger.h"

#include <fstream>
#include <string>

NoiseMSRowProvider::NoiseMap::NoiseMap(std::istream& stream) {
  size_t maxAnt = 0;
  std::string line;
  while (stream) {
    std::getline(stream, line);
    if (stream) {
      std::stringstream linestr(line);
      size_t ant1, ant2;
      linestr >> ant1 >> ant2;
      if (linestr) {
        float stddev;
        linestr >> stddev;
        if (!linestr) stddev = std::numeric_limits<float>::quiet_NaN();
        if (ant1 > ant2) std::swap(ant1, ant2);
        maxAnt = std::max(maxAnt, std::max(ant1, ant2));
        const bool isInserted =
            _map.emplace(std::make_pair(ant1, ant2), stddev).second;
        if (!isInserted)
          throw std::runtime_error(
              "Baseline " + std::to_string(ant1) + " x " +
              std::to_string(ant2) +
              " is specified twice in the noise baseline file");
      }
    }
  }
  Logger::Info << "Read noise baseline file with " << _map.size()
               << " rows and " << maxAnt + 1 << " antennas.\n";
}

float NoiseMSRowProvider::NoiseMap::GetNoiseValue(size_t antenna1,
                                                  size_t antenna2) const {
  size_t a1 = antenna1;
  size_t a2 = antenna2;
  if (a1 > a2) std::swap(a1, a2);
  auto iter = _map.find(std::make_pair(a1, a2));
  if (iter == _map.end())
    throw std::runtime_error(
        "The following baseline was not present in the baseline noise "
        "map: " +
        std::to_string(antenna1) + " x " + std::to_string(antenna2));
  return iter->second;
}

NoiseMSRowProvider::NoiseMSRowProvider(
    const string& msPath, const MSSelection& selection,
    const std::map<size_t, size_t>& selectedDataDescIds,
    const std::string& dataColumnName, bool requireModel)
    : DirectMSRowProvider(msPath, selection, selectedDataDescIds,
                          dataColumnName, requireModel),
      _rng(std::random_device{}()),
      _distribution(0.0, 1.0) {}

void NoiseMSRowProvider::SetNoiseLevel(double noiseStdDevJy) {
  _distribution = std::normal_distribution<float>(0.0, noiseStdDevJy);
}

void NoiseMSRowProvider::SetNoiseBaselineFile(const std::string& filename) {
  std::ifstream file(filename);
  if (!file)
    throw std::runtime_error("Can't open baseline noise file " + filename);
  _noiseMap = NoiseMap(file);
}

void NoiseMSRowProvider::ReadData(DataArray& data, FlagArray& flags,
                                  WeightArray& weights, double& u, double& v,
                                  double& w, uint32_t& dataDescId,
                                  uint32_t& antenna1, uint32_t& antenna2,
                                  uint32_t& fieldId, double& time) {
  DirectMSRowProvider::ReadData(data, flags, weights, u, v, w, dataDescId,
                                antenna1, antenna2, fieldId, time);
  const float stddev =
      _noiseMap.Empty() ? 1.0 : _noiseMap.GetNoiseValue(antenna1, antenna2);

  for (DataArray::contiter iter = data.cbegin(); iter != data.cend(); ++iter) {
    if (std::isfinite(iter->real()) && std::isfinite(iter->imag())) {
      iter->real(_distribution(_rng) * stddev);
      iter->imag(_distribution(_rng) * stddev);
    } else {
      iter->real(std::numeric_limits<float>::quiet_NaN());
      iter->imag(std::numeric_limits<float>::quiet_NaN());
    }
  }
}
