#ifndef PARSET_READER_H
#define PARSET_READER_H

#include <algorithm>
#include <map>
#include <memory>
#include <string>
#include <vector>

#ifdef HAVE_EVERYBEAM
#include <EveryBeam/aterms/parsetprovider.h>
/**
 * @brief Parses the parameter settings (parset) related to the
 * aterm settings in an EveryBeam-acceptable format.
 *
 */
class ParsetReader : public everybeam::aterms::ParsetProvider {
#else
/**
 * @brief Parses the parameter settings (parset) related to the
 * aterm settings.
 *
 */
class ParsetReader {
#endif  // HAVE_EVERYBEAM
 public:
  ParsetReader(const std::string& filename);
  ParsetReader(std::istream& stream);

  class ParsetEntry {
   public:
    ParsetEntry(ParsetEntry&&) = default;
    enum Type { String, StringList };
    ParsetEntry(const std::string& line);
    bool operator<(const ParsetEntry& rhs) const { return _key < rhs._key; }
    const std::string& Key() const { return _key; }
    const std::string& GetStringValue() const;
    const std::vector<std::string>& GetStringListValue() const;

   private:
    struct Value {
      virtual ~Value() {}
    };
    struct StringValue : public Value {
      std::string _value;
    };
    struct StringListValue : public Value {
      std::vector<std::string> _value;
    };

    std::string _key;
    std::unique_ptr<Value> _value;
  };

#ifdef HAVE_EVERYBEAM
  std::string GetString(const std::string& key) const final override;
  std::string GetStringOr(const std::string& key,
                          const std::string& orValue) const final override;
  std::vector<std::string> GetStringList(
      const std::string& key) const final override;

  bool GetBool(const std::string& key) const final override;
  bool GetBoolOr(const std::string& key, bool orValue) const final override;

  double GetDoubleOr(const std::string& key,
                     double orValue) const final override;
#else
  std::string GetString(const std::string& key) const;
  std::string GetStringOr(const std::string& key,
                          const std::string& orValue) const;
  std::vector<std::string> GetStringList(const std::string& key) const;

  bool GetBool(const std::string& key) const;
  bool GetBoolOr(const std::string& key, bool orValue) const;

  double GetDoubleOr(const std::string& key, double orValue) const;
#endif  // HAVE_EVERYBEAM
  // GetDouble is not implemented (and needed) in ParsetProvider
  double GetDouble(const std::string& key) const;

 private:
  void read(std::istream& stream);

  std::map<std::string, ParsetEntry> _entries;
};
#endif
