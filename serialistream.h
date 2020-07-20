#ifndef SERIAL_ISTREAM_H
#define SERIAL_ISTREAM_H

#include <aocommon/uvector.h>
#include "serialostream.h"

#include <stdint.h>

#include <complex>
#include <memory>
#include <vector>

class SerialIStream
{
public:
	typedef std::size_t size_t;
	
	SerialIStream(aocommon::UVector<unsigned char>&& buffer) :
		_buffer(std::move(buffer)),
		_position(_buffer.begin())
	{ }
	
	SerialIStream(SerialOStream&& oStream) :
		_buffer(std::move(oStream._buffer)),
		_position(_buffer.begin())
	{ }
	
	size_t size() const { return _buffer.size(); }
	const unsigned char* data() const { return _buffer.data(); }
	std::string ToString() const { return std::string(reinterpret_cast<const char*>(_buffer.data()), size()); }
	
	const unsigned char* Chunk(size_t size)
	{
		const unsigned char* chunk = &*_position;
		_position += size;
		return chunk;
	}
	
	template<typename T>
	SerialIStream& UInt64(T& value)
	{
		value = (T) read<uint64_t>();
		return *this;
	}
	
	uint64_t UInt64()
	{
		return read<uint64_t>();
	}
	
	template<typename T>
	SerialIStream& UInt32(T& value)
	{
		value = (T) read<uint32_t>();
		return *this;
	}
	
	uint32_t UInt32()
	{
		return read<uint32_t>();
	}
	
	template<typename T>
	SerialIStream& UInt16(T& value)
	{
		value = (T) read<uint16_t>();
		return *this;
	}
	
	uint16_t UInt16()
	{
		return read<uint16_t>();
	}
	
	template<typename T>
	SerialIStream& UInt8(T& value)
	{
		value = (T) read<uint8_t>();
		return *this;
	}
	
	uint8_t UInt8()
	{
		return read<uint8_t>();
	}
	
	SerialIStream& Bool(bool& value)
	{
		value = (UInt8() != 0);
		return *this;
	}
		
	bool Bool()
	{
		return UInt8() != 0;
	}
	
	SerialIStream& Float(float& value)
	{
		return read(value);
	}
	
	float Float()
	{
		return read<float>();
	}
	
	SerialIStream& Double(double& value)
	{
		return read(value);
	}
	
	double Double()
	{
		return read<double>();
	}
	
	SerialIStream& LDouble(long double& value)
	{
		return read(value);
	}
	
	long double LDouble()
	{
		return read<long double>();
	}
	
	SerialIStream& CFloat(std::complex<float>& value)
	{
		return read(value);
	}
	
	std::complex<float> CFloat()
	{
		return read<std::complex<float>>();
	}
	
	SerialIStream& CDouble(std::complex<double>& value)
	{
		return read(value);
	}
	
	std::complex<double> CDouble()
	{
		return read<std::complex<double>>();
	}
	
	SerialIStream& CLDouble(std::complex<long double>& value)
	{
		return read(value);
	}
	
	long double CLDouble()
	{
		return read<long double>();
	}
	
	SerialIStream& String(std::string& str)
	{
		size_t n = UInt64();
		const unsigned char* block = Chunk(n);
		str.resize(n);
		std::copy_n(block, n, str.begin());
		return *this;
	}
	
	std::string String()
	{
		size_t n = UInt64();
		const unsigned char* block = Chunk(n);
		std::string str(n, 0);
		std::copy_n(block, n, str.begin());
		return str;
	}
	
	SerialIStream& VectorUInt64(std::vector<uint64_t>& values)
	{
		// Specialization of VectorUInt64: if T==uint64_t, we don't
		// have to do conversion, so the memory block can be copied directly
		return readVector(values);
	}
	
	template<typename T>
	SerialIStream& VectorUInt64(std::vector<T>& values)
	{
		size_t n = UInt64();
		values.resize(n);
		for(T& val : values)
			UInt64(val);
		return *this;
	}
	
	SerialIStream& VectorFloat(std::vector<float>& values)
	{
		return readVector(values);
	}
	
	SerialIStream& VectorDouble(std::vector<double>& values)
	{
		return readVector(values);
	}
	
	SerialIStream& VectorCFloat(std::vector<std::complex<float>>& values)
	{
		return readVector(values);
	}
	
	SerialIStream& VectorCDouble(std::vector<std::complex<double>>& values)
	{
		return readVector(values);
	}
	
	template<typename T>
	SerialIStream& Ptr(std::unique_ptr<T>& ptr)
	{
		return readPtr(ptr);
	}
	
	template<typename T>
	SerialIStream& Ptr(std::shared_ptr<T>& ptr)
	{
		return readPtr(ptr);
	}
	
	template<typename T>
	SerialIStream& Ptr(T* ptr)
	{
		return readPtr(ptr);
	}
private:
	template<typename T>
	SerialIStream& read(T& value)
	{
		value = *reinterpret_cast<const T*>(Chunk(sizeof(T)));
		return *this;
	}
	
	template<typename T>
	T read()
	{
		return *reinterpret_cast<const T*>(Chunk(sizeof(T)));
	}
	
	template<typename T>
	SerialIStream& readVector(std::vector<T>& values)
	{
		uint64_t size = UInt64();
		values.resize(size);
		size_t n = size * sizeof(T);
		const unsigned char* block = Chunk(n);
		unsigned char* valuePtr = reinterpret_cast<unsigned char*>(values.data());
		std::copy_n(block, n, valuePtr);
		return *this;
	}
	
	template<typename PtrT>
	SerialIStream& readPtr(PtrT& ptr)
	{
		if(Bool())
		{
			ptr.reset(new typename std::remove_reference<decltype(*ptr)>::type());
			ptr->Unserialize(*this);
		}
		else {
			ptr.reset();
		}
		return *this;
	}
	
	aocommon::UVector<unsigned char> _buffer;
	aocommon::UVector<unsigned char>::const_iterator _position;
};

#endif

