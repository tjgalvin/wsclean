#ifndef MODELRENDERER_H
#define MODELRENDERER_H

#include <cstring>

#include <aocommon/polarization.h>

class ModelRenderer {
 public:
  ModelRenderer(long double phaseCentreRA, long double phaseCentreDec,
                long double pixelScaleL, long double pixelScaleM,
                long double phaseCentreDL = 0.0,
                long double phaseCentreDM = 0.0)
      : _phaseCentreRA(phaseCentreRA),
        _phaseCentreDec(phaseCentreDec),
        _pixelScaleL(pixelScaleL),
        _pixelScaleM(pixelScaleM),
        _phaseCentreDL(phaseCentreDL),
        _phaseCentreDM(phaseCentreDM) {}

  /**
   * Restore with circular beam
   */
  void Restore(double* imageData, size_t imageWidth, size_t imageHeight,
               const class Model& model, long double beamSize,
               long double startFrequency, long double endFrequency,
               aocommon::PolarizationEnum polarization);

  /**
   * Restore a model with an elliptical beam
   */
  void Restore(float* imageData, size_t imageWidth, size_t imageHeight,
               const class Model& model, long double beamMaj,
               long double beamMin, long double beamPA,
               long double startFrequency, long double endFrequency,
               aocommon::PolarizationEnum polarization);

  /**
   * Restore elliptical beam using a FFT deconvolution
   */
  void Restore(float* imageData, const float* modelData, size_t imageWidth,
               size_t imageHeight, long double beamMaj, long double beamMin,
               long double beamPA) {
    Restore(imageData, modelData, imageWidth, imageHeight, beamMaj, beamMin,
            beamPA, _pixelScaleL, _pixelScaleM);
  }

  /**
   * Restore elliptical beam using a FFT deconvolution (static version).
   */
  static void Restore(float* imageData, const float* modelData,
                      size_t imageWidth, size_t imageHeight,
                      long double beamMaj, long double beamMin,
                      long double beamPA, long double pixelScaleL,
                      long double pixelScaleM);

  /**
   * Render without beam convolution, such that each point-source is one pixel.
   */
  void RenderModel(float* imageData, size_t imageWidth, size_t imageHeight,
                   const class Model& model, long double startFrequency,
                   long double endFrequency,
                   aocommon::PolarizationEnum polarization) {
    renderModel(imageData, imageWidth, imageHeight, model, startFrequency,
                endFrequency, polarization);
  }

 private:
  void renderPointComponent(float* imageData, size_t imageWidth,
                            size_t imageHeight, long double posRA,
                            long double posDec, long double flux);

  void renderGaussianComponent(float* imageData, size_t imageWidth,
                               size_t imageHeight, long double posRA,
                               long double posDec, long double gausMaj,
                               long double gausMin, long double gausPA,
                               long double flux);

  void renderModel(float* imageData, size_t imageWidth, size_t imageHeight,
                   const class Model& model, long double startFrequency,
                   long double endFrequency,
                   aocommon::PolarizationEnum polarization);

  long double _phaseCentreRA;
  long double _phaseCentreDec;
  long double _pixelScaleL, _pixelScaleM;
  long double _phaseCentreDL, _phaseCentreDM;
  template <typename T>
  static T gaus(T x, T sigma);

  ModelRenderer(const ModelRenderer&) {}
  void operator=(const ModelRenderer&){};
};

#endif
