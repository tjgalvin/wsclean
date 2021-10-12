#include "msweightcolumn.h"

#include "msprovider.h"

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

MsWeightColumn::MsWeightColumn(const casacore::MeasurementSet& ms) {
  if (!MSProvider::OpenWeightSpectrumColumn(ms, weight_column_)) {
    weight_column_.reset(new casacore::ArrayColumn<float>(
        ms, casacore::MS::columnName(casacore::MSMainEnums::WEIGHT)));
    type_ = Type::kScalar;
  }
}

void MsWeightColumn::ReadData(MsRowProviderBase::WeightArray& weights,
                              size_t row,
                              const casacore::IPosition& shape) const {
  switch (type_) {
    case Type::kSpectrum:
      weight_column_->get(row, weights, true);
      break;
    case Type::kScalar: {
      casacore::Array<float> scratch_weights;
      weight_column_->get(row, scratch_weights);
      weights.resize(shape);
      MSProvider::ExpandScalarWeights(scratch_weights, weights);
    } break;
  }
}
