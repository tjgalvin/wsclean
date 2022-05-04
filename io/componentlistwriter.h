#ifndef WSCLEAN_IO_COMPONENT_LIST_WRITER_H_
#define WSCLEAN_IO_COMPONENT_LIST_WRITER_H_

#include <radler/radler.h>
#include <radler/work_table.h>
#include <radler/work_table_entry.h>

#include "../structures/primarybeamimageset.h"

class Settings;

/**
 * @brief Class for extracting the component list from the deconvolution
 * algorithm and writing it to a text file on disk - optionally, components are
 * corrected for the primary beam before writing to disk.
 */
class ComponentListWriter {
 public:
  ComponentListWriter(const Settings& settings,
                      std::unique_ptr<radler::WorkTable> table)
      : settings_(settings), deconvolution_table_(std::move(table)) {}

  /**
   * @brief Save source component list to disk.
   */
  void SaveSourceList(const radler::Radler& deconvolution,
                      long double phase_centre_ra,
                      long double phase_centre_dec) const;

  /**
   * @brief Save primary beam corrected source components to disk.
   */
  void SavePbCorrectedSourceList(const radler::Radler& deconvolution,
                                 long double phase_centre_ra,
                                 long double phase_centre_dec) const;

 private:
  void CorrectChannelForPrimaryBeam(radler::ComponentList& list,
                                    const radler::WorkTableEntry& entry) const;

  PrimaryBeamImageSet LoadAveragePrimaryBeam(size_t image_index) const;

  const Settings& settings_;
  std::unique_ptr<radler::WorkTable> deconvolution_table_;
};

#endif
