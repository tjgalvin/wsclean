#ifndef VISIBILITY_WEIGHTING_MODE_ENUM_H
#define VISIBILITY_WEIGHTING_MODE_ENUM_H

/**
 * This specifies several modi how the visibility weights
 * (normally stored in the SPECTRUM_WEIGHT column)
 * are applied to the data.
 */
enum class VisibilityWeightingMode {
  NormalVisibilityWeighting,
  SquaredVisibilityWeighting,
  UnitVisibilityWeighting
};

#endif