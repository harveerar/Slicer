/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkHaralickTextureFeaturesImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2008-10-16 19:33:45 $
  Version:   $Revision: 1.6 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkHaralickTextureFeaturesImageFilter_h
#define __itkHaralickTextureFeaturesImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkImage.h"
#include "itkNumericTraits.h"
#include "itkVector.h"
//#include "itkScalarImageToTextureFeaturesFilter.h"
#include <vector>

namespace itk
{
//namespace Statistics
//{
/** \class HaralickTextureFeaturesImageFilter
 * \brief Wrapper to calculate the various haralick texture features on the whole image.
 *
 * Computes an image where a given output pixel is computed as the
 * entropy of the input pixels in a neighborhood about the corresponding
 * input pixel. Zero pixels can be excluded by setting a mask with SetMask(). 
 * The mask is assumed to be of the same size as the input image.
 * The output Image is a vector image type (using vector pixels) if we want to output different features.
 * If the output image is set to a scalar image type, then the filter just outputs one of the features in the GetOutput 
 * The user can get various output feature images by specifically calling those GetImages.
 *
 * author: Harini Veeraraghavan
 * email: veerarag at ge dot com
 * date: May 2011
 * This work was developed as part of NAMIC
 *
 * \sa Image
 * \sa Neighborhood
 * \sa NeighborhoodOperator
 * \sa NeighborhoodIterator
 * 
 * \ingroup IntensityImageFilters
 */
// OutputImage is of Type Image< itk::Vector<T, NDimension>, dimensions>

template <class TInputImage, class TOutputImage, class TMaskPixelType=short>
class ITK_EXPORT HaralickTextureFeaturesImageFilter :
    public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
  /** Extract dimension from input and output image. */
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);
  itkStaticConstMacro(OutputImageDimension, unsigned int,
                      TOutputImage::ImageDimension);

  /** Convenient typedefs for simplifying declarations. */
  typedef TInputImage  InputImageType;
  typedef typename InputImageType::Pointer InputImagePointer;
  typedef TOutputImage OutputImageType;

  /** Standard class typedefs. */
  typedef HaralickTextureFeaturesImageFilter                                   Self;
  typedef ImageToImageFilter< InputImageType, OutputImageType> Superclass;
  typedef SmartPointer<Self>                                   Pointer;
  typedef SmartPointer<const Self>                             ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(HaralickTextureFeaturesImageFilter, ImageToImageFilter);
  
  /** Image typedef support. */
  typedef typename InputImageType::PixelType               InputPixelType;
  typedef typename OutputImageType::PixelType              OutputPixelType;
  typedef typename NumericTraits<InputPixelType>::RealType InputRealType;
  
  typedef typename InputImageType::RegionType  InputImageRegionType;
  typedef typename OutputImageType::RegionType OutputImageRegionType;

  typedef typename InputImageType::SizeType InputSizeType;

  typedef Image<TMaskPixelType, itkGetStaticConstMacro(InputImageDimension)>
MaskImageType;
  typedef typename MaskImageType::Pointer MaskImagePointerType;
  typedef TMaskPixelType MaskPixelType;

  //typedef itk::Statistics::ScalarImageToTextureFeaturesFilter< TInputImage, TInputImage > TextureFilterType; 
  
  //itkStaticConstMacro(MeasurementVectorLength, unsigned int, 1);
  //typedef InputPixelType MeasurementType;
  //typedef itk::Vector<MeasurementType, MeasurementVectorLength> MeasurementVectorType;
  //typedef typename itk::Statistics::ListSample<MeasurementVectorType> ListSampleType;
  //typedef itk::Statistics::Histogram< float, itk::Statistics::DenseFrequencyContainer2 > HistogramType;
  //typedef itk::Statistics::SampleToHistogramFilter< ListSampleType, HistogramType > SampleToHistogramFilterType;

  //typedef float HistogramMeasurementType;
  //typedef itk::Statistics::ListSampleToHistogramGenerator<ListSampleType, HistogramMeasurementType, itk::Statistics::DenseFrequencyContainer, MeasurementVectorLength> GeneratorType;

  /** Set the radius of the neighborhood used to compute the Entropy.
   * In general, the radius should be small.  But if set to one, the
   * confidence in the estimates will be marginal. */
  itkSetMacro(Radius, InputSizeType);

  /** Get the radius of the neighborhood used to compute the Entropy */
  itkGetConstReferenceMacro(Radius, InputSizeType);


  itkSetMacro( WindowSize, unsigned long);
  itkGetMacro( WindowSize, unsigned long);
   
  /** Set number of histogram bins and marginal scale for entropy calculation */
  itkSetMacro(NumberOfBins, unsigned long);
  itkGetMacro(NumberOfBins, unsigned long);
  
  /** Set/Get the parameters for normalizing the textures **/
  itkSetMacro( NormalizeTextures, bool);
  itkGetMacro( NormalizeTextures, bool);
  itkBooleanMacro( NormalizeTextures);

  /** Set/Get parameters for median filtering **/
  itkSetMacro( ApplyMedianFiltering, bool);
  itkGetMacro( ApplyMedianFiltering, bool);
  itkBooleanMacro( ApplyMedianFiltering);

  /** Set/Get parameters for normalize range **/
  itkSetMacro( NormalizeMax, float);
  itkGetMacro( NormalizeMax, float);

  itkSetMacro( NormalizeMin, float);
  itkGetMacro( NormalizeMin, float);

  /** Calculate entropy of values stored in listSample */
  //double CalculateEntropy(typename ListSampleType::Pointer listSample);
  
  /** Set the input mask image */
  itkSetObjectMacro(Mask, MaskImageType);
  
  /** Get the input mask image */
  itkGetObjectMacro(Mask, MaskImageType);

   /** output methods for getting the different texture features **/
   const InputImagePointer GetEnergyImage();
  
   const InputImagePointer GetEntropyImage();

   const InputImagePointer GetCorrelationImage();
 
   const InputImagePointer GetDifferenceMomentImage();

   const InputImagePointer GetInertiaImage();

   const InputImagePointer GetClusterShadeImage();

   const InputImagePointer GetClusterProminenceImage();


  /** HaralickTextureFeaturesImageFilter needs a larger input requested region than
   * the output requested region.  As such, HaralickTextureFeaturesImageFilter needs
   * to provide an implementation for GenerateInputRequestedRegion()
   * in order to inform the pipeline execution model.
   *
   * \sa ImageToImageFilter::GenerateInputRequestedRegion() */
  virtual void GenerateInputRequestedRegion() throw(InvalidRequestedRegionError);

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(InputHasNumericTraitsCheck,
                  (Concept::HasNumericTraits<InputPixelType>));
  /** End concept checking */
#endif

protected:
  HaralickTextureFeaturesImageFilter();
  virtual ~HaralickTextureFeaturesImageFilter() {}
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** HaralickTextureFeaturesImageFilter can be implemented as a multithreaded filter.
   * Therefore, this implementation provides a ThreadedGenerateData()
   * routine which is called for each processing thread. The output
   * image data is allocated automatically by the superclass prior to
   * calling ThreadedGenerateData().  ThreadedGenerateData can only
   * write to the portion of the output image specified by the
   * parameter "outputRegionForThread"
   *
   * \sa ImageToImageFilter::ThreadedGenerateData(),
   *     ImageToImageFilter::GenerateData() */
  void ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
                            int threadId );
  //void GenerateData();

  void BeforeThreadedGenerateData();

  void AfterThreadedGenerateData();

  //void Initialize();
private:
  HaralickTextureFeaturesImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  InputSizeType m_Radius;
  unsigned long m_WindowSize;
  unsigned long m_NumberOfBins;
//  unsigned long m_MarginalScale;
  MaskImagePointerType m_Mask;

  bool            m_NormalizeTextures;
  
  bool            m_ApplyMedianFiltering;

  float           m_NormalizeMax;

  float           m_NormalizeMin;

  float  m_MeanEnergy;
  float  m_MeanEntropy;
  float  m_MeanCorrelation;
  float  m_MeanDifferenceMoment;
  float  m_MeanInertia;
  float  m_MeanClusterShade;
  float  m_MeanClusterProminence;

  float  m_StdEnergy;
  float  m_StdEntropy;
  float  m_StdCorrelation;
  float  m_StdDifferenceMoment;
  float  m_StdInertia;
  float  m_StdClusterShade;
  float  m_StdClusterProminence;

  float  m_KurtosisEnergy;
  float  m_KurtosisEntropy;
  float  m_KurtosisCorrelation;
  float  m_KurtosisDifferenceMoment;
  float  m_KurtosisInertia;
  float  m_KurtosisClusterShade;
  float  m_KurtosisClusterProminence;


   

};
  
} // end namespace itk
//} // end namespace Statistics

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkHaralickTextureFeaturesImageFilter.txx"
#endif

#endif
