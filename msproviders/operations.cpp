#include "operations.h"

#include <aocommon/logger.h>

#include <casacore/tables/DataMan/DataManager.h>
#include <casacore/tables/Tables/ArrColDesc.h>

using aocommon::Logger;
using schaapcommon::reordering::MSSelection;
using schaapcommon::reordering::StorageManagerType;

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

void GetRowRange(casacore::MeasurementSet& ms, const MSSelection& selection,
                 size_t& startRow, size_t& endRow) {
  startRow = 0;
  endRow = ms.nrow();
  if (selection.HasInterval()) {
    Logger::Info << "Determining first and last row index... ";
    Logger::Info.Flush();
    casacore::ScalarColumn<double> timeColumn(
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

bool OpenWeightSpectrumColumn(
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

void InitializeModelColumn(casacore::MeasurementSet& ms,
                           const std::string& data_column_name,
                           const std::string& model_column_name,
                           StorageManagerType type) {
  casacore::ArrayColumn<casacore::Complex> data_column(ms, data_column_name);
  ms.reopenRW();
  if (ms.tableDesc().isColumn(model_column_name)) {
    casacore::ArrayColumn<casacore::Complex> model_column(ms,
                                                          model_column_name);
    const bool is_defined = model_column.isDefined(0);
    bool is_same_shape = false;
    if (is_defined) {
      if (type == StorageManagerType::Sisco) {
        is_same_shape = true;
      } else {
        casacore::IPosition model_shape = model_column.shape(0);
        casacore::IPosition data_shape = data_column.shape(0);
        is_same_shape = model_shape.isEqual(data_shape);
      }
    }
    if (!is_defined || !is_same_shape) {
      Logger::Warn << "WARNING: Your model column does not have the same shape "
                      "as your data column: resetting MODEL column.\n";
      FillModelColumn(data_column, model_column);
    }
  } else {  // No column exists with the given model_column_name
    Logger::Info << "Adding model data column " << model_column_name << "... ";
    Logger::Info.Flush();
    std::string st_man_name = "StandardStMan";
    bool use_direct_column = false;
    casacore::Record dataman_specification;
    switch (type) {
      case StorageManagerType::Default:
        break;
      case StorageManagerType::StokesI:
        st_man_name = "StokesIStMan";
        use_direct_column = true;
        break;
      case StorageManagerType::Sisco:
        st_man_name = "SiscoStMan";
        dataman_specification.define("predict_level", 2);
        dataman_specification.define("deflate_level", 9);
        break;
    }
    casacore::DataManagerCtor constructor =
        casacore::DataManager::getCtor(st_man_name);
    std::unique_ptr<casacore::DataManager> st_man(
        constructor(model_column_name + "_dm", dataman_specification));
    if (!st_man)
      throw std::runtime_error(
          st_man_name +
          " storage manager requested, but it is not available in "
          "casacore");
    casacore::ArrayColumnDesc<casacore::Complex> model_column_desc(
        model_column_name);
    if (use_direct_column) {
      model_column_desc.setShape(data_column.shape(0));
      model_column_desc.setOptions(casacore::ColumnDesc::Direct |
                                   casacore::ColumnDesc::FixedShape);
    }
    casacore::TableDesc table_desc;
    table_desc.addColumn(model_column_desc, model_column_name);
    ms.addColumn(table_desc, *st_man, true);

    casacore::ArrayColumn<casacore::Complex> model_column(ms,
                                                          model_column_name);
    if (type != StorageManagerType::Sisco)
      FillModelColumn(data_column, model_column);

    Logger::Info << "DONE\n";
  }
}

casacore::ArrayColumn<float> InitializeImagingWeightColumn(
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

std::set<aocommon::PolarizationEnum> GetMSPolarizations(
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

double GetMsInterval(const casacore::MeasurementSet& ms) {
  const casacore::ScalarColumn<double> interval_column(
      ms, ms.columnName(casacore::MSMainEnums::INTERVAL));
  return interval_column(0);
}

std::vector<std::string> GetAntennaNames(const casacore::MSAntenna& antenna) {
  const casacore::ScalarColumn<casacore::String> antennaNameColumn(
      antenna, antenna.columnName(casacore::MSAntenna::NAME));

  std::vector<std::string> antenna_names;
  antenna_names.reserve(antennaNameColumn.nrow());
  for (size_t i = 0; i < antennaNameColumn.nrow(); ++i) {
    antenna_names.push_back(antennaNameColumn(i));
  }
  return antenna_names;
}

}  // namespace wsclean
