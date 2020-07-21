/***************************************************************************
 *   Copyright (C) 2008 by A.R. Offringa   *
 *   offringa@astro.rug.nl   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#ifndef SERIALIZABLE_H
#define SERIALIZABLE_H

#include <aocommon/uvector.h>

#include <stdint.h>

#include <iostream>
#include <complex>
#include <vector>

class Serializable {
 public:
  virtual ~Serializable() {}
  virtual void Serialize(std::ostream& stream) const = 0;
  virtual void Unserialize(std::istream& stream) = 0;

  template <typename T>
  static void SerializeToUInt64(std::ostream& stream, T value) {
    uint64_t val64t = value;
    stream.write(reinterpret_cast<char*>(&val64t), sizeof(val64t));
  }

  template <typename T>
  static void SerializeToUInt32(std::ostream& stream, T value) {
    uint32_t val32t = value;
    stream.write(reinterpret_cast<char*>(&val32t), sizeof(val32t));
  }

  template <typename T>
  static void SerializeToUInt16(std::ostream& stream, T value) {
    uint16_t val16t = value;
    stream.write(reinterpret_cast<char*>(&val16t), sizeof(val16t));
  }

  template <typename T>
  static void SerializeToUInt8(std::ostream& stream, T value) {
    uint8_t val8t = value;
    stream.write(reinterpret_cast<char*>(&val8t), sizeof(val8t));
  }

  static void SerializeToFloat(std::ostream& stream, float value) {
    stream.write(reinterpret_cast<char*>(&value), sizeof(value));
  }

  static void SerializeToFloatC(std::ostream& stream,
                                std::complex<float> value) {
    stream.write(reinterpret_cast<char*>(&value), sizeof(value));
  }

  static void SerializeToDouble(std::ostream& stream, double value) {
    stream.write(reinterpret_cast<char*>(&value), sizeof(value));
  }

  static void SerializeToDoubleC(std::ostream& stream,
                                 std::complex<double> value) {
    stream.write(reinterpret_cast<char*>(&value), sizeof(value));
  }

  static void SerializeToLDouble(std::ostream& stream, long double value) {
    stream.write(reinterpret_cast<char*>(&value), sizeof(value));
  }

  static void SerializeToLDoubleC(std::ostream& stream,
                                  std::complex<long double> value) {
    stream.write(reinterpret_cast<char*>(&value), sizeof(value));
  }

  static void SerializeToBool(std::ostream& stream, bool value) {
    SerializeToUInt8(stream, value ? 1 : 0);
  }

  static void SerializeToString(std::ostream& stream, const std::string& str) {
    SerializeToUInt64(stream, str.size());
    stream.write(str.c_str(), str.size());
  }

  static uint64_t UnserializeUInt64(std::istream& stream) {
    return Unserialize<uint64_t>(stream);
  }

  static uint32_t UnserializeUInt32(std::istream& stream) {
    return Unserialize<uint32_t>(stream);
  }

  static uint16_t UnserializeUInt16(std::istream& stream) {
    return Unserialize<uint16_t>(stream);
  }

  static uint8_t UnserializeUInt8(std::istream& stream) {
    return Unserialize<uint8_t>(stream);
  }

  static bool UnserializeBool(std::istream& stream) {
    return UnserializeUInt8(stream) != 0;
  }

  static float UnserializeFloat(std::istream& stream) {
    return Unserialize<float>(stream);
  }

  static std::complex<float> UnserializeFloatC(std::istream& stream) {
    return Unserialize<std::complex<float>>(stream);
  }

  static double UnserializeDouble(std::istream& stream) {
    return Unserialize<double>(stream);
  }

  static std::complex<double> UnserializeDoubleC(std::istream& stream) {
    return Unserialize<std::complex<double>>(stream);
  }

  static long double UnserializeLDouble(std::istream& stream) {
    return Unserialize<long double>(stream);
  }

  static std::complex<long double> UnserializeLDoubleC(std::istream& stream) {
    return Unserialize<std::complex<long double>>(stream);
  }

  static void UnserializeString(std::istream& stream, std::string& destStr) {
    size_t size = UnserializeUInt64(stream);
    aocommon::UVector<char> str(size);
    stream.read(str.data(), size);
    destStr = std::string(str.data(), size);
  }

  static std::string UnserializeString(std::istream& stream) {
    size_t size = UnserializeUInt64(stream);
    aocommon::UVector<char> str(size);
    stream.read(str.data(), size);
    return std::string(str.data(), size);
  }

  template <typename T>
  static void SerializeVectorUInt64(std::ostream& stream,
                                    const std::vector<T>& values) {
    uint64_t size = values.size();
    SerializeToUInt64(stream, size);
    for (const T val : values) SerializeToUInt64(stream, val);
  }

  template <typename T>
  static void UnserializeVectorUInt64(std::istream& stream,
                                      std::vector<T>& values) {
    uint64_t size = UnserializeUInt64(stream);
    values.resize(size);
    for (size_t i = 0; i != size; ++i) values[i] = UnserializeUInt64(stream);
  }

  static void SerializeVectorFloat(std::ostream& stream,
                                   const std::vector<float>& values) {
    SerializeToUInt64(stream, values.size());
    for (const float& val : values) SerializeToFloat(stream, val);
  }

  static std::vector<float> UnserializeVectorFloat(std::istream& stream) {
    std::vector<float> values(UnserializeUInt64(stream));
    for (float& val : values) val = UnserializeFloat(stream);
    return values;
  }

  static void SerializeVectorDouble(std::ostream& stream,
                                    const std::vector<double>& values) {
    SerializeToUInt64(stream, values.size());
    for (const double& val : values) SerializeToDouble(stream, val);
  }

  static std::vector<double> UnserializeVectorDouble(std::istream& stream) {
    std::vector<double> values(UnserializeUInt64(stream));
    for (double& val : values) val = UnserializeDouble(stream);
    return values;
  }

  static void SerializeVectorFloatC(
      std::ostream& stream, const std::vector<std::complex<float>>& values) {
    SerializeToUInt64(stream, values.size());
    for (const std::complex<float>& val : values)
      SerializeToFloatC(stream, val);
  }

  static std::vector<std::complex<float>> UnserializeVectorFloatC(
      std::istream& stream) {
    std::vector<std::complex<float>> values(UnserializeUInt64(stream));
    for (std::complex<float>& val : values) val = UnserializeFloatC(stream);
    return values;
  }

  template <typename T>
  static void SerializePtr(std::ostream& stream,
                           const std::unique_ptr<T>& ptr) {
    if (ptr) {
      Serializable::SerializeToBool(stream, true);
      ptr->Serialize(stream);
    } else {
      Serializable::SerializeToBool(stream, false);
    }
  }

  template <typename T>
  static void UnserializePtr(std::istream& stream, std::unique_ptr<T>& ptr) {
    bool hasValue = Serializable::UnserializeBool(stream);
    if (hasValue) {
      ptr.reset(new T());
      ptr->Unserialize(stream);
    } else {
      ptr.reset();
    }
  }

  template <typename T>
  static void SerializePtr(std::ostream& stream,
                           const std::shared_ptr<T>& ptr) {
    if (ptr) {
      Serializable::SerializeToBool(stream, true);
      ptr->Serialize(stream);
    } else {
      Serializable::SerializeToBool(stream, false);
    }
  }

  template <typename T>
  static void UnserializePtr(std::istream& stream, std::shared_ptr<T>& ptr) {
    bool hasValue = Serializable::UnserializeBool(stream);
    if (hasValue) {
      ptr.reset(new T());
      ptr->Unserialize(stream);
    } else {
      ptr.reset();
    }
  }

 private:
  template <typename T>
  static T Unserialize(std::istream& stream) {
    T val;
    stream.read(reinterpret_cast<char*>(&val), sizeof(val));
    return val;
  }
};

#endif
