#ifndef MSPROVIDER_H
#define MSPROVIDER_H

#include "synchronizedms.h"

#include "../structures/observationinfo.h"

#include <aocommon/multibanddata.h>
#include <aocommon/polarization.h>

#include <complex>
#include <optional>
#include <set>
#include <type_traits>
#include <vector>

namespace casacore {
class MeasurementSet;
}  // namespace casacore

namespace schaapcommon::reordering {
struct ChannelRange;
class MSSelection;
}  // namespace schaapcommon::reordering

namespace wsclean {
class MSReader;

aocommon::MultiBandData MakeSelectedPartBands(
    const aocommon::MultiBandData& input,
    const aocommon::VectorMap<schaapcommon::reordering::ChannelRange>& ranges);

std::vector<aocommon::MultiBandData> MakeSelectedBands(
    const aocommon::MultiBandData& input,
    const std::vector<
        aocommon::VectorMap<schaapcommon::reordering::ChannelRange>>& channels);

bool HasFrequencyBda(const casacore::MeasurementSet& ms);

/**
 * The abstract MSProvider class is the base class for classes that read and
 * write the visibilities. Write functionality is directly provided by classes
 * that derive from MSProvider, whereas reading functionality is provided in a
 * separate class (MSReader), which can be instantiated with @ref
 * MSProvider::MakeReader. An MSProvider knows which rows are selected and
 * doesn't write to unselected rows. This information on selected rows is also
 * passed to the @ref MSReader. MSProvider provides the visibilities weighted
 * with the visibility weight and converts the visibilities to a requested
 * polarization. The @ref ContiguousMS and
 * @ref ReorderedMs classes implement the MSProvider interface.
 *
 * The class maintains an index for the write position. The index for the
 * reading position is maintained by the closely connected @ref MSReader class.
 * Writing (and reading) goes sequentially through the data.
 *
 * An MS provider provides data for a single dataDescID, and therefore for
 * a single spectral window. Measurement sets with multiple data desc IDs
 * require creating multiple MS providers. Because of this, all rows of a
 * single MS provider have the same number of channels.
 */
class MSProvider {
 public:
  struct MetaData {
    double u_in_m;
    double v_in_m;
    double w_in_m;
    double time;
    uint32_t data_desc_id;
    uint32_t field_id;
    uint32_t antenna1;
    uint32_t antenna2;
  };

  MSProvider() = default;

  virtual ~MSProvider();

  MSProvider(const MSProvider&) = delete;
  MSProvider& operator=(const MSProvider&) = delete;

  /**
   * Returns a descriptive string that identifies the part. In the case of a
   * on-disk ms provider, this could e.g. be the filename + part index. It
   * is used by WSClean to notify the user what it is working on.
   */
  virtual std::string PartDescription() const = 0;

  /**
   * Move the model writing position to the next row.
   */
  virtual void NextOutputRow() = 0;

  /**
   * Reset the writing position to the first row.
   */
  virtual void ResetWritePosition() = 0;

  /**
   * Write model visibilities to the current writing position. If add is true,
   * the provided data are add-assigned to the existing model visibilities. If
   * false, the existing model visibilities are overwritten.
   */
  virtual void WriteModel(const std::complex<float>* buffer, bool add) = 0;

  /**
   * Prepare the msprovider for writing. This is explicitly required before
   * writing to make it possible to image read-only sets when no writing is
   * required.
   */
  virtual void ReopenRW() = 0;

  /**
   * First MJD time of the data.
   */
  virtual double StartTime() = 0;

  /**
   * Polarization that this msprovider provides.
   * Be aware that it may return 'DiagonalInstrumental' or 'Instrumental',
   * which means it provides 2 or 4 polarizations that are provided by the
   * underlying measurement set / data.
   */
  virtual aocommon::PolarizationEnum Polarization() = 0;

  /**
   * Number of channels provided by this provider. If the set is regular,
   * this is equal to the number of channels in every row. This value may be
   * different from the underlying measurement set if not all channels are
   * selected.
   */
  virtual size_t NMaxChannels() = 0;

  /**
   * Does every row have the same number of channels? If true, fixed size
   * arrays can be allocated beforehand for reading/writing. This is also
   * used in reordering to determine whether a size-per-row needs to be
   * stored.
   *
   * Both BDA or the use of multiple spectral windows with different number of
   * channels can make the data irregular. Even if a MS has regular subbands,
   * the data can become irregular because of selection or output channel
   * splitting.
   *
   * In the reordering, a stricter definition is used: being regular
   * means there that a provider is strictly for one spectral window. This means
   * that a provider for two SPWs with the same nr of channels would still
   * return 'false'.
   */
  virtual bool IsRegular() const = 0;

  /**
   * Are data_desc_ids in this measurement sets part of a BDA averaging
   * scheme? If true, @ref IsRegular() will also always return false. If this
   * methods returns false, the ms may still have multiple data_desc_ids, but
   * these then represent genuine bands instead of different averaging factors.
   * In this situation, @ref IsRegular() may still return true if those bands
   * have different number of channels.
   */
  virtual bool HasFrequencyBda() const = 0;

  /**
   * Count of antennas in the underlying measurement set (irrespective of
   * whether all of them are selected). An antenna indix returned by the @ref
   * ReadMeta() methods is always less than this number.
   */
  virtual size_t NAntennas() = 0;

  /**
   * This is the number of polarizations provided by this MSProvider.
   * Note that this does not have to be equal to the nr of pol in the
   * MS, as in most cases each pol is provided by a separate msprovider.
   */
  virtual size_t NPolarizations() = 0;

  virtual size_t NRows() = 0;

  /**
   * Integration time of one timestep, in seconds. If time BDA is used, then it
   * is undefined what is returned.
   */
  virtual double Interval() = 0;

  virtual ObservationInfo GetObservationInfo() = 0;

  /**
   * When applying facet solutions, the antenna names of the msprovider should
   * match the names given in h5parm solution files.
   */
  virtual std::vector<std::string> GetAntennaNames() = 0;

  /**
   * Get the band information that this MSProvider covers after applying
   * any selection. The indexing of channels matches therefore with
   * the indexing of visibilities and weights by reading/writing.
   */
  virtual const aocommon::MultiBandData& SelectedBands() = 0;

  virtual std::unique_ptr<MSReader> MakeReader() = 0;

  /**
   * Returns the underlying measurement set if it is available. Requiring an
   * on-disk measurement set should be avoided as much as possible, since not
   * all providers have an on-disk measurement set (e.g. the InMemoryProvider).
   * Nevertheless, certain very specific features may require access to the
   * MS to read in non-standard data. As of yet, applying the beam is one of
   * those cases.
   */
  virtual std::optional<SynchronizedMS> MsIfAvailable() { return {}; }

  /**
   * Reset model data in the MSProvider to zeros.
   */
  virtual void ResetModelColumn();

  static void CopyRealToComplex(std::complex<float>* dest, const float* source,
                                size_t n) {
    const float* end = source + n;
    while (source != end) {
      *dest = *source;
      ++dest;
      ++source;
    }
  }
};

}  // namespace wsclean

#endif
