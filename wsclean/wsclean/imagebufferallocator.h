#ifndef IMAGE_BUFFER_ALLOCATOR_H
#define IMAGE_BUFFER_ALLOCATOR_H

#include <complex>
#include <iostream>
#include <vector>
#include <stdexcept>
#include <mutex>
#include "logger.h"

//#define USE_DIRECT_ALLOCATOR

#ifndef USE_DIRECT_ALLOCATOR

/**
 * The ImageBufferAllocator is a simple allocator for maintaining large memory chunks.
 * 
 * This class was made to avoid segmentation: it was found that wsclean initially used more
 * memory than it had allocated, which was due to iteratively freeing and allocating
 * large chunks, but because of intermediate small allocations the large chunks no longer
 * fitted, and new memory was returned.
 */
class ImageBufferAllocator
{
public:
	/** A unique smart pointer to an image buffer. */
	template<typename NumT>
	class TPtr
	{
	public:
		friend class ImageBufferAllocator;
		/** Construct empty pointer */
		TPtr()
			: _data(nullptr), _allocator(nullptr) { }
		TPtr(std::nullptr_t)
			: _data(nullptr), _allocator(nullptr) { }
		/** Construct initialized pointer
			* @param data databuffer
			* @param allocator corresponding allocator
			*/
		TPtr(NumT* data, ImageBufferAllocator& allocator) 
			: _data(data), _allocator(&allocator) { }
			
		TPtr(TPtr<NumT>&& source)
		  : _data(source._data), _allocator(source._allocator)
		{
			source._data = nullptr;
			source._allocator = nullptr;
		}
		/** Destructor */
		~TPtr()
		{
			reset();
		}
		TPtr<NumT>& operator=(TPtr<NumT>&& rhs)
		{
			reset(rhs._data, *rhs._allocator);
			rhs._data = nullptr;
			rhs._data = nullptr;
			return *this;
		}
		bool operator==(std::nullptr_t) const { return _data == nullptr; }
		/** Is the ptr set? */
		operator bool() const { return _data != nullptr; }
		/** @returns pointer to data */
		NumT* data() const { return _data; }
		/** @returns reference to data */ 
		NumT& operator*() const { return *_data; }
		/** Destructs buffer if set. */
		void reset()
		{
			// We have to check for nullptr, because the _allocator might not be set.
			if(_data != nullptr) _allocator->Free(_data);
			_data = nullptr;
		}
		/** Destructs buffer if set and assign new contents.
			* @param data databuffer
			* @param allocator corresponding allocator
			*/
		void reset(NumT* data, ImageBufferAllocator& allocator)
		{
			// We have to check for nullptr, because the _allocator might not be set.
			if(_data != nullptr) _allocator->Free(_data);
			_data = data;
			_allocator = &allocator;
		}
		/** Retrieve image element with given index. */
		NumT& operator[](size_t index) const { return _data[index]; }
		
	private:
		TPtr(const TPtr<NumT>&) = delete;
		void operator=(const TPtr<NumT>&) = delete;
		
		NumT* _data;
		ImageBufferAllocator* _allocator;
	};
	
	typedef TPtr<double> Ptr;
	
	/** A unique smart pointer to a complex image buffer. */
	template<typename NumT>
	class CPtr
	{
	public:
		friend class ImageBufferAllocator;
		/** Construct empty pointer */
		CPtr()
			: _data(nullptr), _allocator(nullptr) { }
		CPtr(std::nullptr_t)
			: _data(nullptr), _allocator(nullptr) { }
		/** Construct initialized pointer
			* @param data databuffer
			* @param allocator corresponding allocator
			*/
		CPtr(std::complex<NumT>* data, ImageBufferAllocator& allocator) 
			: _data(data), _allocator(&allocator) { }
			
		CPtr(CPtr<NumT>&& source)
		  : _data(source._data), _allocator(source._allocator)
		{
			source._data = nullptr;
			source._allocator = nullptr;
		}
		/** Destructor */
		~CPtr()
		{
			reset();
		}
		CPtr<NumT>& operator=(CPtr<NumT>&& rhs)
		{
			reset(rhs._data, *rhs._allocator);
			rhs._data = nullptr;
			rhs._data = nullptr;
			return *this;
		}
		bool operator==(std::nullptr_t) const { return _data == nullptr; }
		/** Is the ptr set? */
		operator bool() const { return _data != nullptr; }
		/** @returns pointer to data */
		std::complex<NumT>* data() const { return _data; }
		/** @returns reference to data */ 
		std::complex<NumT>& operator*() const { return *_data; }
		/** Destructs buffer if set. */
		void reset()
		{
			// We have to check for nullptr, because the _allocator might not be set.
			if(_data != nullptr) _allocator->Free(_data);
			_data = nullptr;
		}
		/** Destructs buffer if set and assign new contents.
			* @param data databuffer
			* @param allocator corresponding allocator
			*/
		void reset(std::complex<NumT>* data, ImageBufferAllocator& allocator)
		{
			// We have to check for nullptr, because the _allocator might not be set.
			if(_data != nullptr) _allocator->Free(_data);
			_data = data;
			_allocator = &allocator;
		}
		/** Retrieve image element with given index. */
		std::complex<NumT>& operator[](size_t index) const { return _data[index]; }
		
	private:
		CPtr(const CPtr<NumT>&) = delete;
		void operator=(const CPtr<NumT>&) = delete;
		
		std::complex<NumT>* _data;
		ImageBufferAllocator* _allocator;
	};
	
	ImageBufferAllocator() : _buffers(), _nReal(0), _nComplex(0), _nRealMax(0), _nComplexMax(0), _previousSize(0), _direct(false)
	{ }
	
	ImageBufferAllocator(bool direct) : _buffers(), _nReal(0), _nComplex(0), _nRealMax(0), _nComplexMax(0), _previousSize(0), _direct(direct)
	{ }
	
	ImageBufferAllocator(const ImageBufferAllocator& source) = delete;
	
	ImageBufferAllocator& operator=(const ImageBufferAllocator& source) = delete;
	
	ImageBufferAllocator& operator=(ImageBufferAllocator&& source)
	{
		clear();
		std::swap(_buffers, source._buffers);
		std::swap(_nReal, source._nReal);
		std::swap(_nComplex, source._nComplex);
		std::swap(_nRealMax, source._nRealMax);
		std::swap(_nComplexMax, source._nComplexMax);
		std::swap(_previousSize, source._previousSize);
		std::swap(_direct, source._direct);
		return *this;
	}
	
	~ImageBufferAllocator()
	{
		clear();
	}
	
	void ReportStatistics() const
	{
		std::lock_guard<std::mutex> guard(_mutex);
		double totalSize = 0.0;
		for(const Buffer& buffer : _buffers)
		{
			totalSize += double(buffer.size) * double(sizeof(double)*2);
		}
		Logger::Info << "Image buf alloc stats:\n"
			"         max alloc'd images = " << _nRealMax << " real + " << _nComplexMax << " complex\n"
			"       max allocated chunks = " << _buffers.size() << "\n"
			"      current allocated mem = " << round(totalSize/1e8)/10.0 << " GB \n";
	}
	
	void Allocate(size_t size, Ptr& ptr)
	{
		ptr.reset(Allocate(size), *this);
	}
	
	Ptr AllocatePtr(size_t size)
	{
		return Ptr(Allocate(size), *this);
	}
	
	template<typename NumT>
	TPtr<NumT> AllocateTPtr(size_t size);
	
	double* Allocate(size_t size)
	{
		if(_direct)
			return directAllocate(size);
		else {
			if(size == 0) size=2;
			else if(size%2==1) ++size;
			std::lock_guard<std::mutex> guard(_mutex);
			
			if(size != _previousSize)
			{
				FreeUnused();
				_previousSize = size;
			}
			
			++_nReal;
			if(_nReal > _nRealMax) _nRealMax = _nReal;
			for(std::vector<Buffer>::iterator i=_buffers.begin(); i!=_buffers.end(); ++i)
			{
				if(i->size == size)
				{
					if(!i->isFirstHalfUsed)
					{
						i->isFirstHalfUsed = true;
						return i->ptr;
					}
					if(!i->isSecondHalfUsed)
					{
						i->isSecondHalfUsed = true;
						return i->ptr + size;
					}
				}
			}
			Buffer* newBuffer = allocateNewBuffer(size);
			newBuffer->isFirstHalfUsed = true;
			return newBuffer->ptr;
		}
	}
	
	template<typename NumT>
	CPtr<NumT> AllocateCPtr(size_t size);
	
	std::complex<double>* AllocateComplex(size_t size)
	{
		if(_direct)
			return reinterpret_cast<std::complex<double>*>(directAllocate(size*2));
		else {
			if(size == 0) size=2;
			else if(size%2==1) ++size;
			std::lock_guard<std::mutex> guard(_mutex);
			
			if(size != _previousSize)
			{
				FreeUnused();
				_previousSize = size;
			}
			
			++_nComplex;
			if(_nComplex > _nComplexMax) _nComplexMax = _nComplex;
			for(std::vector<Buffer>::iterator i=_buffers.begin(); i!=_buffers.end(); ++i)
			{
				if(i->size == size)
				{
					if(!i->isFirstHalfUsed && !i->isSecondHalfUsed)
					{
						i->isFirstHalfUsed = true;
						i->isSecondHalfUsed = true;
						return reinterpret_cast<std::complex<double>*>(i->ptr);
					}
				}
			}
			Buffer* newBuffer = allocateNewBuffer(size);
			newBuffer->isFirstHalfUsed = true;
			newBuffer->isSecondHalfUsed = true;
			return reinterpret_cast<std::complex<double>*>(newBuffer->ptr);
		}
	}
	
	void Free(float* buffer)
	{
		Free(reinterpret_cast<double*>(buffer));
	}
	
	void Free(double* buffer)
	{
		if(_direct)
			free(buffer);
		else if(buffer != nullptr)
		{
			std::lock_guard<std::mutex> guard(_mutex);
			bool found = false;
			for(std::vector<Buffer>::iterator i=_buffers.begin(); i!=_buffers.end(); ++i)
			{
				if(i->ptr == buffer)
				{
					found = true;
					i->isFirstHalfUsed = false;
					break;
				}
				else if(i->ptr + i->size == buffer)
				{
					found = true;
					i->isSecondHalfUsed = false;
					break;
				}
			}
			if(!found)
			{
				std::cerr << "Invalid or double call to ImageBufferAllocator::Free(double*).\n";
				throw std::runtime_error("Invalid or double call to ImageBufferAllocator::Free(double*).");
			}
			--_nReal;
		}
	}
	
	void Free(std::complex<float>* buffer)
	{ Free(reinterpret_cast<std::complex<double>*>(buffer)); }
	
	void Free(std::complex<double>* buffer)
	{
		if(_direct)
			free(buffer);
		else if(buffer != nullptr)
		{
			std::lock_guard<std::mutex> guard(_mutex);
			bool found = false;
			for(std::vector<Buffer>::iterator i=_buffers.begin(); i!=_buffers.end(); ++i)
			{
				if(i->ptr == reinterpret_cast<double*>(buffer))
				{
					found = true;
					i->isFirstHalfUsed = false;
					i->isSecondHalfUsed = false;
					break;
				}
			}
			if(!found)
			{
				std::cerr << "Invalid or double call to ImageBufferAllocator::Free(std::complex<double>*).\n";
				throw std::runtime_error("Invalid or double call to ImageBufferAllocator::Free(std::complex<double>*).");
			}
			--_nComplex;
		}
	}
	
	void FreeUnused()
	{
		size_t unusedCount = 0;
		std::vector<Buffer>::iterator i=_buffers.begin();
		while(i!=_buffers.end())
		{
			if(!i->isFirstHalfUsed && !i->isSecondHalfUsed)
			{
				free(i->ptr);
				_buffers.erase(i);
				i = _buffers.begin();
				++unusedCount;
			} else {
				++i;
			}
		}
		if(unusedCount != 0)
		{
			Logger::Debug << "Freed " << unusedCount << " image buffer(s).\n";
		}
	}
	
private:
	struct Buffer
	{
		double* ptr;
		size_t size;
		bool isFirstHalfUsed, isSecondHalfUsed;
	};
	
	Buffer* allocateNewBuffer(size_t size)
	{
		_buffers.push_back(Buffer());
		Buffer* buffer = &_buffers.back();
		buffer->ptr = directAllocate(size*2);
		buffer->size = size;
		buffer->isFirstHalfUsed = false;
		buffer->isSecondHalfUsed = false;
		return buffer;
	}
	
	double* directAllocate(size_t size)
	{
		double* ptr;
		int errVal = posix_memalign(reinterpret_cast<void**>(&ptr), sizeof(double), size*sizeof(double));
		if(errVal != 0)
		{
			std::ostringstream msg;
			msg << "posix_memalign() failed when allocating " << size*sizeof(double) * 2 << " bytes: ";
			switch(errVal)
			{
			case EINVAL:
				msg << "the alignment argument was not a power of two, or was not a multiple of sizeof(void *)";
				break;
			case ENOMEM:
				msg << "there was insufficient memory to fulfill the allocation request.";
				break;
			default:
				msg << "an unknown error value was returned";
				break;
			}
			throw std::runtime_error(msg.str());
		}
		return ptr;
	}
	
	void clear()
	{
		std::lock_guard<std::mutex> guard(_mutex);
		size_t usedCount = 0;
		std::ostringstream str;
		for(Buffer& buffer : _buffers)
		{
			if(buffer.isFirstHalfUsed) {
				++usedCount;
				str << "Still used: buffer of " << buffer.size << '\n';
				buffer.isFirstHalfUsed = false;
			}
			if(buffer.isSecondHalfUsed) {
				++usedCount;
				str << "Still used: buffer of " << buffer.size << '\n';
				buffer.isSecondHalfUsed = false;
			}
			// In case of leak, leave this allocated so that e.g. valgrind can diagnose it as well.
			if(!buffer.isFirstHalfUsed && !buffer.isSecondHalfUsed)
				free(buffer.ptr);
		}
		if(usedCount != 0)
		{
			std::cerr << usedCount << " image buffer(s) were still in use when image buffer allocator was destroyed!\n" << str.str();
		}
	}
	
	std::vector<Buffer> _buffers;
	size_t _nReal, _nComplex, _nRealMax, _nComplexMax, _previousSize;
	bool _direct;
	mutable std::mutex _mutex;
};

template<>
inline ImageBufferAllocator::TPtr<double> ImageBufferAllocator::AllocateTPtr(size_t size)
{
	return AllocatePtr(size);
}

template<>
inline ImageBufferAllocator::TPtr<float> ImageBufferAllocator::AllocateTPtr(size_t size)
{
	return TPtr<float>(reinterpret_cast<float*>(Allocate((size+1)/2)), *this);
}

template<>
inline ImageBufferAllocator::CPtr<double> ImageBufferAllocator::AllocateCPtr(size_t size)
{
	return CPtr<double>(AllocateComplex(size), *this);
}

template<>
inline ImageBufferAllocator::CPtr<float> ImageBufferAllocator::AllocateCPtr(size_t size)
{
	return CPtr<float>(reinterpret_cast<std::complex<float>*>(AllocateComplex((size+1)/2)), *this);
}

#else // USE_DIRECT_ALLOCATOR

class ImageBufferAllocator
{
public:
	class Ptr
	{
	public:
		friend class ImageBufferAllocator;
		/** Construct empty pointer */
		Ptr()
			: _data(nullptr), _allocator(nullptr) { }
		/** Construct initialized pointer
			* @param data databuffer
			* @param allocator corresponding allocator
			*/
		Ptr(double* data, ImageBufferAllocator& allocator) 
			: _data(data), _allocator(&allocator) { }
			
		Ptr(Ptr&& source)
		  : _data(source._data), _allocator(source._allocator)
		{
			source._data = nullptr;
			source._allocator = nullptr;
		}
		/** Destructor */
		~Ptr()
		{
			reset();
		}
		Ptr& operator=(Ptr&& rhs)
		{
			reset(rhs._data, *rhs._allocator);
			rhs._data = nullptr;
			rhs._data = nullptr;
			return *this;
		}
		bool operator==(std::nullptr_t) const { return _data == nullptr; }
		/** Is the ptr set? */
		operator bool() const { return _data != nullptr; }
		/** @returns pointer to data */
		double* data() const { return _data; }
		/** @returns reference to data */ 
		double& operator*() const { return *_data; }
		/** Destructs buffer if set. */
		void reset()
		{
			// We have to check for nullptr, because the _allocator might not be set.
			if(_data != nullptr) _allocator->Free(_data);
			_data = nullptr;
		}
		/** Destructs buffer if set and assign new contents.
			* @param data databuffer
			* @param allocator corresponding allocator
			*/
		void reset(double* data, ImageBufferAllocator& allocator)
		{
			// We have to check for nullptr, because the _allocator might not be set.
			if(_data != nullptr) _allocator->Free(_data);
			_data = data;
			_allocator = &allocator;
		}
		/** Retrieve image element with given index. */
		double& operator[](size_t index) const { return _data[index]; }
	private:
		Ptr(const Ptr&) = delete;
		void operator=(const Ptr&) = delete;
		
		double* _data;
		ImageBufferAllocator* _allocator;
	};
	
	void ReportStatistics() const
	{
	}
	
	double* Allocate(size_t size)
	{
		return new double[size];
	}
	
	void Allocate(size_t size, Ptr& ptr)
	{
		ptr.reset(Allocate(size), *this);
	}
	
	std::complex<double>* AllocateComplex(size_t size)
	{
		return new std::complex<double>[size];
	}
	
	void Free(double* buffer)
	{
		delete[] buffer;
	}
	
	void Free(std::complex<double>* buffer)
	{
		delete[] buffer;
	}
	
	void FreeUnused() { }
};

#endif // USE_DIRECT_ALLOCATOR

#endif
