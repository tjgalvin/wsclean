#ifndef NUMBERLIST_H
#define NUMBERLIST_H

#include <concepts>
#include <cstdlib>
#include <ranges>
#include <string>

#include <aocommon/uvector.h>

namespace wsclean {

template <typename T>
concept NumericRangeLike =
    std::ranges::input_range<T> && requires(std::ranges::range_value_t<T> v) {
      { std::to_string(v) } -> std::convertible_to<std::string>;
    };

class NumberList {
 public:
  static aocommon::UVector<int> ParseIntList(const std::string& str) {
    aocommon::UVector<int> list;
    std::string temp = str;
    size_t pos = temp.find(",");
    while (pos != std::string::npos) {
      std::string idStr = temp.substr(0, pos);
      temp = temp.substr(pos + 1);
      int num = atoi(idStr.c_str());
      list.push_back(num);
      pos = temp.find(",");
    }
    int num = atoi(temp.c_str());
    list.push_back(num);
    return list;
  }

  static aocommon::UVector<double> ParseDoubleList(const std::string& str) {
    aocommon::UVector<double> list;
    std::string temp = str;
    size_t pos = temp.find(",");
    while (pos != std::string::npos) {
      std::string idStr = temp.substr(0, pos);
      temp = temp.substr(pos + 1);
      double num = atof(idStr.c_str());
      list.push_back(num);
      pos = temp.find(",");
    }
    double num = atof(temp.c_str());
    list.push_back(num);
    return list;
  }

  template <NumericRangeLike T>
  static std::string ToString(const T& data) {
    if (std::empty(data)) return "";
    std::string result;
    if constexpr (std::ranges::sized_range<T>) {
      result.reserve(std::ranges::size(data) * 5);
    }
    bool first = true;
    for (const auto& element : data) {
      if (!first) result += ", ";
      result += std::to_string(element);
      first = false;
    }
    return result;
  }
};

}  // namespace wsclean

#endif
