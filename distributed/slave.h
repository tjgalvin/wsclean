#ifndef SLAVE_H
#define SLAVE_H

#include "../wsclean/wscleansettings.h"

class Slave {
 public:
  Slave(const WSCleanSettings& settings) : _settings(settings) {}

  void Run();

 private:
  void grid(size_t bodySize);

  const WSCleanSettings _settings;
};

#endif
