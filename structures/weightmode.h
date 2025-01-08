#ifndef STRUCTURES_WEIGHTMODE_H_
#define STRUCTURES_WEIGHTMODE_H_

#include <aocommon/io/serialostream.h>
#include <aocommon/io/serialistream.h>

#include <string>
#include <sstream>

namespace wsclean {

enum class WeightClass { Natural, Uniform, Briggs };

class WeightMode {
 public:
  explicit constexpr WeightMode(WeightClass weight_class)
      : class_(weight_class) {}

  constexpr WeightMode(const WeightMode& source)
      : class_(source.class_),
        briggs_robustness_(source.briggs_robustness_),
        super_weight_(source.super_weight_) {}

  WeightMode& operator=(const WeightMode& source) = default;

  constexpr bool operator==(const WeightMode& rhs) const {
    if (class_ != rhs.class_ || super_weight_ != rhs.super_weight_)
      return false;
    else if (class_ == WeightClass::Briggs)
      return briggs_robustness_ == rhs.briggs_robustness_;
    else
      return true;
  }

  constexpr static WeightMode Briggs(double briggsRobustness) {
    WeightMode m(WeightClass::Briggs);
    m.briggs_robustness_ = briggsRobustness;
    return m;
  }

  constexpr WeightClass Class() const { return class_; }
  constexpr bool IsNatural() const { return class_ == WeightClass::Natural; }
  constexpr bool IsUniform() const { return class_ == WeightClass::Uniform; }
  constexpr bool IsBriggs() const { return class_ == WeightClass::Briggs; }

  constexpr double BriggsRobustness() const { return briggs_robustness_; }
  constexpr double SuperWeight() const { return super_weight_; }

  constexpr void SetSuperWeight(double superWeight) {
    super_weight_ = superWeight;
  }

  constexpr bool RequiresGridding() const {
    return true;
  }  // { IsUniform() || IsBriggs(); }

  std::string ToString() const {
    switch (class_) {
      case WeightClass::Uniform:
        return "uniform";
      case WeightClass::Natural:
        return "natural";
      case WeightClass::Briggs: {
        std::ostringstream s;
        s << "Briggs'(" << briggs_robustness_ << ")";
        return s.str();
      }
      default:
        return "?";
    }
  }

  void Serialize(aocommon::SerialOStream& stream) const {
    stream.UInt32(static_cast<uint32_t>(class_));
    if (class_ == WeightClass::Briggs) {
      stream.Double(briggs_robustness_);
      stream.Double(super_weight_);
    }
  }

  void Unserialize(aocommon::SerialIStream& stream) {
    class_ = WeightClass(stream.UInt32());
    if (class_ == WeightClass::Briggs) {
      stream.Double(briggs_robustness_);
      stream.Double(super_weight_);
    }
  }

 private:
  WeightClass class_;
  double briggs_robustness_ = 0.0;
  double super_weight_ = 1.0;
};

}  // namespace wsclean

#endif
