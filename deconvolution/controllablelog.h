#ifndef CONTROLLABLE_LOG_H
#define CONTROLLABLE_LOG_H

#include "../wsclean/logger.h"

#include <mutex>
#include <string>
#include <vector>

class ControllableLog : public LogReceiver {
 public:
  ControllableLog(std::mutex* mutex)
      : _mutex(mutex), _atNewLine(true), _isMuted(false), _isActive(true) {}

  virtual ~ControllableLog() {}

  void Mute(bool mute) { _isMuted = mute; }
  bool IsMuted() const { return _isMuted; }

  void Activate(bool active) { _isActive = active; }
  bool IsActive() const { return _isActive; }

  void SetTag(const std::string& tag) { _tag = tag; }
  void SetOutputOnce(const std::string& str) { _outputOnce = str; }

 protected:
  virtual void output(enum Logger::LoggerLevel level,
                      const std::string& str) final override {
    if (!str.empty()) {
      std::lock_guard<std::mutex> lock(*_mutex);

      bool skip =
          ((level == Logger::DebugLevel || level == Logger::InfoLevel) &&
           _isMuted) ||
          (level == Logger::DebugLevel && !Logger::IsVerbose());

      if (!skip) {
        _lineBuffer += str;
        if (_lineBuffer.back() == '\n') {
          if (!_outputOnce.empty()) {
            forward(level, _outputOnce);
            _outputOnce.clear();
          }
          forward(level, _tag);
          forward(level, _lineBuffer);
          _lineBuffer.clear();
        }
      }
    }
  }

 private:
  std::mutex* _mutex;
  std::string _tag;
  bool _atNewLine, _isMuted, _isActive;
  std::string _lineBuffer, _outputOnce;
};

class FacetLogSet {
 public:
  FacetLogSet() : _nHorizontal(0), _nVertical(0) {}

  void Initialize(size_t nHorizontal, size_t nVertical) {
    std::lock_guard<std::mutex> lock(_outputMutex);
    _nHorizontal = nHorizontal;
    _nVertical = nVertical;
    size_t n = nHorizontal * nVertical;
    _logs.clear();
    _logs.reserve(n);
    for (size_t i = 0; i != n; ++i) {
      _logs.emplace_back(&_outputMutex);
      _logs[i].SetTag("P" + std::to_string(i) + ": ");
      _logs[i].Mute(true);
      _logs[i].Activate(false);
    }
  }

  ControllableLog& operator[](size_t index) { return _logs[index]; }

  void Activate(size_t index) {
    std::lock_guard<std::mutex> oLock(_outputMutex);
    if (!_logs[index].IsActive()) {
      _logs[index].Activate(true);

      unmuteMostCentral();
    }
  }

  void Deactivate(size_t index) {
    std::lock_guard<std::mutex> oLock(_outputMutex);
    if (_logs[index].IsActive()) {
      _logs[index].Mute(true);
      _logs[index].SetOutputOnce(std::string());
      _logs[index].Activate(false);

      unmuteMostCentral();
    }
  }

 private:
  void unmuteMostCentral() {
    size_t unmutedLog = _logs.size();
    for (size_t i = 0; i != _logs.size(); ++i) {
      if (!_logs[i].IsMuted()) unmutedLog = i;
      _logs[i].SetTag("P" + std::to_string(i) + ": ");
      _logs[i].Mute(true);
    }

    // Find an active facet that is as close as possible to
    // the centre (since these often contain the most interesting info)
    bool found = false;
    size_t fX = 0, fY = 0, fD = 0;
    for (size_t y = 0; y != _nVertical; ++y) {
      for (size_t x = 0; x != _nHorizontal; ++x) {
        size_t index = y * _nHorizontal + x;
        if (_logs[index].IsActive()) {
          int dx = int(x) - int(_nHorizontal) / 2,
              dy = int(y) - int(_nVertical) / 2;
          size_t distSq = dx * dx + dy * dy;
          if (!found || distSq < fD) {
            fX = x;
            fY = y;
            fD = distSq;
            found = true;
          }
        }
      }
    }

    if (found) {
      size_t index = fY * _nHorizontal + fX;
      _logs[index].Mute(false);
      _logs[index].SetTag(" ");

      if (index != unmutedLog) {
        _logs[index].SetOutputOnce("Switching to output of subimage " +
                                   std::to_string(index) + "\n");
      }
    }
  }

  std::mutex _outputMutex;
  std::vector<ControllableLog> _logs;
  size_t _nHorizontal, _nVertical;
};

#endif
