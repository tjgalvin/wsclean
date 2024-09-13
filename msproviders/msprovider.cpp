#include "msprovider.h"
#include "msreaders/msreader.h"

#include <aocommon/logger.h>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/tables/Tables/ScalarColumn.h>
#include <casacore/tables/Tables/ArrColDesc.h>

#include "../structures/msselection.h"

using aocommon::Logger;

namespace wsclean {
namespace {

void FillModelColumn(const casacore::ArrayColumn<casacore::Complex>& dataColumn,
                     casacore::ArrayColumn<casacore::Complex>& modelColumn) {
  casacore::Array<casacore::Complex> zeroArray;
  for (size_t row = 0; row != dataColumn.nrow(); ++row) {
    zeroArray.resize(dataColumn.shape(row));
    zeroArray = casacore::Complex(0.0, 0.0);
    modelColumn.put(row, zeroArray);
  }
}
}  // namespace

MSProvider::~MSProvider() {}

void MSProvider::GetRowRange(casacore::MeasurementSet& ms,
                             const MSSelection& selection, size_t& startRow,
                             size_t& endRow) {
  startRow = 0;
  endRow = ms.nrow();
  if (selection.HasInterval()) {
    Logger::Info << "Determining first and last row index... ";
    Logger::Info.Flush();
    casacore::ROScalarColumn<double> timeColumn(
        ms, casacore::MS::columnName(casacore::MSMainEnums::TIME));
    double time = timeColumn(0);
    size_t timestepIndex = 0;
    for (size_t row = 0; row != ms.nrow(); ++row) {
      if (time != timeColumn(row)) {
        ++timestepIndex;
        if (timestepIndex == selection.IntervalStart()) startRow = row;
        if (timestepIndex == selection.IntervalEnd()) {
          endRow = row;
          break;
        }
        time = timeColumn(row);
      }
    }
    Logger::Info << "DONE (" << startRow << '-' << endRow << ")\n";
  }
}

void MSProvider::GetRowRangeAndIDMap(casacore::MeasurementSet& ms,
                                     const MSSelection& selection,
                                     size_t& startRow, size_t& endRow,
                                     const std::set<size_t>& dataDescIds,
                                     std::vector<size_t>& idToMSRow) {
  startRow = 0;
  endRow = ms.nrow();

  Logger::Info << "Mapping measurement set rows... ";
  Logger::Info.Flush();
  casacore::ArrayColumn<double> uvwColumn(
      ms, casacore::MS::columnName(casacore::MSMainEnums::UVW));
  casacore::ScalarColumn<int> antenna1Column(
      ms, casacore::MS::columnName(casacore::MSMainEnums::ANTENNA1));
  casacore::ScalarColumn<int> antenna2Column(
      ms, casacore::MS::columnName(casacore::MSMainEnums::ANTENNA2));
  casacore::ScalarColumn<int> fieldIdColumn(
      ms, casacore::MS::columnName(casacore::MSMainEnums::FIELD_ID));
  casacore::ScalarColumn<double> timeColumn(
      ms, casacore::MS::columnName(casacore::MSMainEnums::TIME));
  casacore::ScalarColumn<int> dataDescIdColumn(
      ms, ms.columnName(casacore::MSMainEnums::DATA_DESC_ID));
  double time = timeColumn(0);
  size_t timestepIndex = 0;
  bool timeStepSelected =
      !selection.HasInterval() || timestepIndex == selection.IntervalStart();
  for (size_t row = 0; row != ms.nrow(); ++row) {
    if (time != timeColumn(row)) {
      ++timestepIndex;
      if (selection.HasInterval() &&
          timestepIndex == selection.IntervalStart()) {
        startRow = row;
        timeStepSelected = true;
      }
      if (timestepIndex == selection.IntervalEnd()) {
        if (selection.HasInterval()) endRow = row;
        break;
      }
      time = timeColumn(row);
    }
    if (timeStepSelected) {
      const int a1 = antenna1Column(row), a2 = antenna2Column(row),
                fieldId = fieldIdColumn(row),
                dataDescId = dataDescIdColumn(row);
      casacore::Vector<double> uvw = uvwColumn(row);
      std::set<size_t>::const_iterator dataDescIdIter =
          dataDescIds.find(dataDescId);
      if (selection.IsSelected(fieldId, timestepIndex, a1, a2, uvw.data()) &&
          dataDescIdIter != dataDescIds.end())
        idToMSRow.push_back(row);
    }
  }
  Logger::Info << "DONE (" << startRow << '-' << endRow << "; "
               << idToMSRow.size() << " rows)\n";
}

void MSProvider::InitializeModelColumn(casacore::MeasurementSet& ms) {
  casacore::ArrayColumn<casacore::Complex> dataColumn(
      ms, casacore::MS::columnName(casacore::MSMainEnums::DATA));
  ms.reopenRW();
  if (ms.isColumn(casacore::MSMainEnums::MODEL_DATA)) {
    casacore::ArrayColumn<casacore::Complex> modelColumn(
        ms, casacore::MS::columnName(casacore::MSMainEnums::MODEL_DATA));
    bool isDefined = modelColumn.isDefined(0);
    bool isSameShape = false;
    if (isDefined) {
      casacore::IPosition modelShape = modelColumn.shape(0);
      casacore::IPosition dataShape = dataColumn.shape(0);
      isSameShape = modelShape == dataShape;
    }
    if (!isDefined || !isSameShape) {
      Logger::Warn << "WARNING: Your model column does not have the same shape "
                      "as your data column: resetting MODEL column.\n";
      FillModelColumn(dataColumn, modelColumn);
    }
  } else {  // No column named MODEL_DATA
    Logger::Info << "Adding model data column... ";
    Logger::Info.Flush();
    casacore::ArrayColumnDesc<casacore::Complex> modelColumnDesc(
        ms.columnName(casacore::MSMainEnums::MODEL_DATA));
    try {
      ms.addColumn(modelColumnDesc, "StandardStMan", true, true);
    } catch (std::exception& e) {
      ms.addColumn(modelColumnDesc, "StandardStMan", false, true);
    }

    casacore::ArrayColumn<casacore::Complex> modelColumn(
        ms, casacore::MS::columnName(casacore::MSMainEnums::MODEL_DATA));
    FillModelColumn(dataColumn, modelColumn);

    Logger::Info << "DONE\n";
  }
}

casacore::ArrayColumn<float> MSProvider::InitializeImagingWeightColumn(
    casacore::MeasurementSet& ms) {
  ms.reopenRW();
  casacore::ArrayColumn<casacore::Complex> dataColumn(
      ms, casacore::MS::columnName(casacore::MSMainEnums::DATA));
  if (ms.tableDesc().isColumn("IMAGING_WEIGHT_SPECTRUM")) {
    return casacore::ArrayColumn<float>(ms, "IMAGING_WEIGHT_SPECTRUM");
  } else {
    Logger::Info << "Adding imaging weight spectrum column... ";
    Logger::Info.Flush();
    casacore::IPosition shape = dataColumn.shape(0);
    casacore::ArrayColumnDesc<float> modelColumnDesc("IMAGING_WEIGHT_SPECTRUM",
                                                     shape);
    try {
      ms.addColumn(modelColumnDesc, "StandardStMan", true, true);
    } catch (std::exception& e) {
      ms.addColumn(modelColumnDesc, "StandardStMan", false, true);
    }

    casacore::Array<float> zeroArray(shape);
    for (casacore::Array<float>::contiter i = zeroArray.cbegin();
         i != zeroArray.cend(); ++i)
      *i = 0.0;

    casacore::ArrayColumn<float> imgWColumn(ms, "IMAGING_WEIGHT_SPECTRUM");
    for (size_t row = 0; row != ms.nrow(); ++row)
      imgWColumn.put(row, zeroArray);
    Logger::Info << "DONE\n";
    return imgWColumn;
  }
}

std::set<aocommon::PolarizationEnum> MSProvider::GetMSPolarizations(
    size_t data_desc_id, const casacore::MeasurementSet& ms) {
  // First get the polarization index corresponding with the data desc id
  casacore::MSDataDescription data_description_table = ms.dataDescription();
  casacore::ScalarColumn<int> polarization_index_column(
      data_description_table,
      casacore::MSDataDescription::columnName(
          casacore::MSDataDescription::POLARIZATION_ID));
  const size_t polarization_index = polarization_index_column(data_desc_id);
  casacore::MSPolarization pol_table = ms.polarization();
  std::set<aocommon::PolarizationEnum> pols;
  casacore::ArrayColumn<int> corr_type_column(
      pol_table, casacore::MSPolarization::columnName(
                     casacore::MSPolarizationEnums::CORR_TYPE));

  // Now get the information corresponding with the polarization index
  casacore::Array<int> corr_type_vec(corr_type_column(polarization_index));
  for (casacore::Array<int>::const_contiter p = corr_type_vec.cbegin();
       p != corr_type_vec.cend(); ++p) {
    pols.emplace(aocommon::Polarization::AipsIndexToEnum(*p));
  }

  return pols;
}

void MSProvider::ResetModelColumn() {
  std::unique_ptr<MSReader> msReader = MakeReader();
  const std::vector<std::complex<float>> buffer(NChannels() * NPolarizations(),
                                                {0.0f, 0.0f});
  while (msReader->CurrentRowAvailable()) {
    // Always overwrite
    const bool addToMS = false;
    WriteModel(buffer.data(), addToMS);
    NextOutputRow();
    msReader->NextInputRow();
  }
}

bool MSProvider::OpenWeightSpectrumColumn(
    const casacore::MeasurementSet& ms,
    std::unique_ptr<casacore::ArrayColumn<float>>& weightColumn) {
  bool isWeightDefined;
  if (ms.isColumn(casacore::MSMainEnums::WEIGHT_SPECTRUM)) {
    weightColumn.reset(new casacore::ArrayColumn<float>(
        ms, casacore::MS::columnName(casacore::MSMainEnums::WEIGHT_SPECTRUM)));
    isWeightDefined = weightColumn->isDefined(0);
  } else {
    isWeightDefined = false;
  }
  if (!isWeightDefined) {
    Logger::Warn
        << "WARNING: This measurement set has no or an invalid WEIGHT_SPECTRUM "
           "column; will use less informative WEIGHT column.\n";
    weightColumn.reset();
  }
  return isWeightDefined;
}

}  // namespace wsclean
