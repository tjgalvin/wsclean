#ifndef CONTROLLABLE_LOG_H
#define CONTROLLABLE_LOG_H

#include <aocommon/logger.h>

#include <mutex>
#include <string>
#include <vector>

class ControllableLog final : public aocommon::LogReceiver {
 public:
  ControllableLog(std::mutex* mutex)
      : _mutex(mutex), _isMuted(false), _isActive(true) {}

  virtual ~ControllableLog() {}

  void Mute(bool mute) { _isMuted = mute; }
  bool IsMuted() const { return _isMuted; }

  void Activate(bool active) { _isActive = active; }
  bool IsActive() const { return _isActive; }

  void SetTag(const std::string& tag) { _tag = tag; }
  void SetOutputOnce(const std::string& str) { _outputOnce = str; }

 private:
  void Output(enum aocommon::Logger::LoggerLevel level,
              const std::string& str) override {
    if (!str.empty()) {
      std::lock_guard<std::mutex> lock(*_mutex);

      bool skip = ((level == aocommon::Logger::kDebugLevel ||
                    level == aocommon::Logger::kInfoLevel) &&
                   _isMuted) ||
                  (level == aocommon::Logger::kDebugLevel &&
                   !aocommon::Logger::IsVerbose());

      if (!skip) {
        _lineBuffer += str;
        if (_lineBuffer.back() == '\n') {
          if (!_outputOnce.empty()) {
            Forward(level, _outputOnce);
            _outputOnce.clear();
          }
          Forward(level, _tag);
          Forward(level, _lineBuffer);
          _lineBuffer.clear();
        }
      }
    }
  }

  std::mutex* _mutex;
  std::string _tag;
  bool _isMuted;
  bool _isActive;
  std::string _lineBuffer;
  std::string _outputOnce;
};

#endif
