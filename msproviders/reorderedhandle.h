#ifndef WSCLEAN_REORDERED_HANDLE_H_
#define WSCLEAN_REORDERED_HANDLE_H_

#include <memory>

#include <aocommon/io/serialistream.h>
#include <aocommon/io/serialostream.h>

#include <schaapcommon/reordering/handledata.h>

namespace wsclean {

class ReorderedHandle {
  // Both these are a friend of ReorderedHandle
  // in order to access the data_ member.
  friend class ReorderedMsReader;
  friend class ReorderedMsProvider;

 public:
  ReorderedHandle() = default;

  ReorderedHandle(std::unique_ptr<schaapcommon::reordering::HandleData> data)
      : data_(std::move(data)) {}

  void Serialize(aocommon::SerialOStream& stream) const { stream.Ptr(data_); }

  void Unserialize(aocommon::SerialIStream& stream) { stream.Ptr(data_); }

 private:
  std::shared_ptr<schaapcommon::reordering::HandleData> data_;
};

}  // namespace wsclean

#endif
