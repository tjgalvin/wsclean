#ifndef MULTIBANDDATA_H
#define MULTIBANDDATA_H

#include <stdexcept>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/Tables/ScalarColumn.h>

#include <aocommon/banddata.h>

using aocommon::BandData;

/**
 * Contains information about a set of bands. This follows the CASA Measurement
 * Set model; one MultiBandData instance can contain the band information
 * contained in the CASA Measurement Set.
 */
class MultiBandData {
 public:
  using iterator = std::vector<BandData>::iterator;
  using const_iterator = std::vector<BandData>::const_iterator;

  /**
   * Construct an empty MultiBandData.
   */
  MultiBandData() {}

  /**
   * Construct a MultiBandData from a Measurement Set.
   * @param ms A measurement set. MultiBandData reads the spectral window table
   * and the data description table of this measurement set.
   */
  explicit MultiBandData(const casacore::MeasurementSet& ms)
      : MultiBandData(ms.spectralWindow(), ms.dataDescription()) {}

  /**
   * Construct a MultiBandData from the Measurement Set tables.
   * @param spwTable The spectral window table of a measurement set.
   * @param dataDescTable The data description table of a measurement set.
   */
  MultiBandData(const casacore::MSSpectralWindow& spwTable,
                const casacore::MSDataDescription& dataDescTable);

  /**
   * Construct a MultiBandData from another instance but only select a part of
   * each band data. This function also works when not all bands have the
   * same number of channels. If endChannel is larger than the number of
   * channels for one of the bands, the band is selected up to its last channel.
   * @param source Other instance that will be partially copied.
   * @param startChannel Start of channel range to initialize this instance
   * with.
   * @param endChannel End of channel range (exclusive) to initialize this
   * instance with.
   */
  MultiBandData(const MultiBandData& source, size_t startChannel,
                size_t endChannel);

  /**
   * Index operator to retrieve a band data given a dataDescID.
   * @param dataDescID A valid data description ID for which the band is
   * returned.
   * @returns The BandData for the requested band.
   */
  const BandData& operator[](size_t dataDescID) const {
    return _bandData[_dataDescToBand[dataDescID]];
  }

  /**
   * Get number of bands stored.
   * @returns Number of bands.
   */
  size_t BandCount() const { return _bandData.size(); }

  /**
   * Returns the unique number of data description IDs.
   * @returns Unique number of data desc IDs.
   */
  size_t DataDescCount() const { return _dataDescToBand.size(); }

  /**
   * Get lowest frequency.
   * @returns The channel frequency of the channel with lowest frequency.
   */
  double LowestFrequency() const {
    if (_bandData.empty()) return 0.0;
    double freq = _bandData[0].LowestFrequency();
    for (size_t i = 0; i != _bandData.size(); ++i)
      freq = std::min(freq, _bandData[i].LowestFrequency());
    return freq;
  }

  /**
   * Get centre frequency.
   * @returns (BandStart() + BandEnd()) * 0.5.
   */
  double CentreFrequency() const { return (BandStart() + BandEnd()) * 0.5; }

  /**
   * Get highest frequency.
   * @returns The channel frequency of the channel with highest frequency.
   */
  double HighestFrequency() const {
    if (_bandData.empty()) return 0.0;
    double freq = _bandData[0].HighestFrequency();
    for (size_t i = 0; i != _bandData.size(); ++i)
      freq = std::max(freq, _bandData[i].HighestFrequency());
    return freq;
  }

  /**
   * Get total bandwidth covered.
   * @returns BandEnd() - BandStart().
   */
  double Bandwidth() const { return BandEnd() - BandStart(); }

  /**
   * Get the start frequency of the lowest frequency channel.
   * @return Start of covered bandwidth.
   */
  double BandStart() const {
    if (_bandData.empty()) return 0.0;
    double freq = std::min(_bandData[0].BandStart(), _bandData[0].BandEnd());
    for (size_t i = 0; i != _bandData.size(); ++i)
      freq = std::min(
          freq, std::min(_bandData[i].BandStart(), _bandData[i].BandEnd()));
    return freq;
  }

  /**
   * Get the end frequency of the highest frequency channel.
   * @return End of covered bandwidth.
   */
  double BandEnd() const {
    if (_bandData.empty()) return 0.0;
    double freq = std::max(_bandData[0].BandStart(), _bandData[0].BandEnd());
    for (size_t i = 0; i != _bandData.size(); ++i)
      freq = std::max(
          freq, std::max(_bandData[i].BandStart(), _bandData[i].BandEnd()));
    return freq;
  }

  /**
   * Get the maximum number of channels in a band.
   * @returns Maximum number of channels.
   */
  size_t MaxChannels() const {
    size_t maxChannels = 0;
    for (const BandData& band : _bandData) {
      if (band.ChannelCount() > maxChannels) maxChannels = band.ChannelCount();
    }
    return maxChannels;
  }

  /**
   * Map a dataDescId to the corresponding band index.
   * @param dataDescId A dataDescId as e.g. used in a main table.
   * @returns The band index, which is equal to the row index in the spw
   * table that describes the band in a measurement set.
   */
  size_t GetBandIndex(size_t dataDescId) const {
    return _dataDescToBand[dataDescId];
  }

  /**
   * Compose a list of dataDescIds that are used in the measurement set.
   * "Used" here means it is references from the main table.
   * @param mainTable the measurement set.
   * @returns Set of used dataDescIds.
   */
  std::set<size_t> GetUsedDataDescIds(
      casacore::MeasurementSet& mainTable) const;

  /**
   * Adds a new band at the end of the list of bands.
   * The band will be linked to the first available dataDescId, which
   * is the number returned by @ref DataDescCount().
   * @returns the dataDescId of this band.
   */
  size_t AddBand(const BandData& data) {
    const size_t swdId = _dataDescToBand.size();
    const size_t bandId = _bandData.size();
    _dataDescToBand.emplace_back(bandId);
    _bandData.emplace_back(data);
    return swdId;
  }

  iterator begin() { return _bandData.begin(); }
  const_iterator begin() const { return _bandData.begin(); }

  iterator end() { return _bandData.end(); }
  const_iterator end() const { return _bandData.end(); }

 private:
  std::vector<size_t> _dataDescToBand;
  std::vector<BandData> _bandData;
};

#endif
