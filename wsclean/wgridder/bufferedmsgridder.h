#ifndef BUFFERED_MS_GRIDDER_H
#define BUFFERED_MS_GRIDDER_H

#include "../wsclean/msgridderbase.h"

#include "../lane.h"
#include "../multibanddata.h"

#include <complex>
#include <memory>

#include <casacore/casa/Arrays/Array.h>
#include <casacore/tables/Tables/ArrayColumn.h>

#include <thread>

namespace casacore {
	class MeasurementSet;
}
class ImageBufferAllocator;

class BufferedMSGridder : public MSGridderBase
{
	public:
		BufferedMSGridder(class ImageBufferAllocator* imageAllocator, size_t threadCount, double memFraction, double absMemLimit);
	
		virtual void Invert() final override;
		
		virtual void Predict(ImageBufferAllocator::Ptr image) final override;
		virtual void Predict(ImageBufferAllocator::Ptr, ImageBufferAllocator::Ptr) final override
		{
			throw std::runtime_error("Can not do imaginary imaging in this mode");
		}
		
		virtual ImageBufferAllocator::Ptr ImageRealResult() final override
		{
			return std::move(_image);
		}
		virtual ImageBufferAllocator::Ptr ImageImaginaryResult() final override
		{
			throw std::runtime_error("Can not do imaginary imaging in this mode");
		}
		
		virtual bool HasGriddingCorrectionImage() const final override { return false; }
		virtual void GetGriddingCorrectionImage(double *) const final override { }
		
		virtual size_t ActualInversionWidth() const final override { return _actualInversionWidth; }
		virtual size_t ActualInversionHeight() const final override { return _actualInversionHeight; }
		
		virtual void FreeImagingData() final override
		{ }
		
		virtual size_t getSuggestedWGridSize() const final override { return 1; }
		
	private:
		ImageBufferAllocator::Ptr _image;
		
		void gridMeasurementSet(MSData& msData);

		void predictMeasurementSet(MSData& msData);

		size_t _cpuCount;
		int64_t _memSize;
		ImageBufferAllocator* _imageBufferAllocator;
		std::unique_ptr<class WGriddingGridder_Simple> _gridder;
};

#endif
