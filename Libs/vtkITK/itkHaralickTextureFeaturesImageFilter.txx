/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkHaralickTextureFeaturesImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2009-04-06 00:19:17 $
  Version:   $Revision: 1.7 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkHaralickTextureFeaturesImageFilter_txx
#define __itkHaralickTextureFeaturesImageFilter_txx

#include "itkHaralickTextureFeaturesImageFilter.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkNeighborhoodIterator.h"
#include "itkImageRegionIterator.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkConstantBoundaryCondition.h"
#include "itkOffset.h"
#include "itkProgressReporter.h"
#include "itkImageFileWriter.h"

#include "itkNumericTraits.h"
#include "itkVectorContainer.h"
#include "itkVector.h"

#include "itkRegionOfInterestImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkScalarImageToTextureFeaturesFilter.h"

#include "itkMedianImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"

#include "itkVectorIndexSelectionCastImageFilter.h"

#include <ctime>
#include <vcl_cmath.h>
#include <vcl_algorithm.h>

namespace itk
{

template <class TInputImage, class TOutputImage, class TMaskPixelType>
HaralickTextureFeaturesImageFilter<TInputImage, TOutputImage, TMaskPixelType>
::HaralickTextureFeaturesImageFilter()
{
  m_Radius.Fill(1);
  m_NumberOfBins = 128; //64;
  m_Mask = NULL;
  m_WindowSize = 3;

  // normalize parameters
  m_NormalizeTextures = true;
  m_ApplyMedianFiltering = false;

  m_NormalizeMax = 255.0;
  m_NormalizeMin = 0.0;

}

template <class TInputImage, class TOutputImage, class TMaskPixelType>
void 
HaralickTextureFeaturesImageFilter<TInputImage, TOutputImage, TMaskPixelType>
::GenerateInputRequestedRegion() throw (InvalidRequestedRegionError)
{
  // call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();
  
  // get pointers to the input and output
  typename Superclass::InputImagePointer inputPtr = 
    const_cast< TInputImage * >( this->GetInput() );
  typename Superclass::OutputImagePointer outputPtr = this->GetOutput();
  
  if ( !inputPtr || !outputPtr )
    {
    return;
    }

  // get a copy of the input requested region (should equal the output
  // requested region)
  typename TInputImage::RegionType inputRequestedRegion;
  inputRequestedRegion = inputPtr->GetRequestedRegion();

  // pad the input requested region by the operator radius
  inputRequestedRegion.PadByRadius( m_Radius );

  // crop the input requested region at the input's largest possible region
  if ( inputRequestedRegion.Crop(inputPtr->GetLargestPossibleRegion()) )
    {
    inputPtr->SetRequestedRegion( inputRequestedRegion );
    inputPtr->SetBufferedRegion( inputRequestedRegion );
    return;
    }
  else
    {
    // Couldn't crop the region (requested region is outside the largest
    // possible region).  Throw an exception.

    // store what we tried to request (prior to trying to crop)
    inputPtr->SetRequestedRegion( inputRequestedRegion );
    inputPtr->SetBufferedRegion( inputRequestedRegion );
    
    // build an exception
    InvalidRequestedRegionError e(__FILE__, __LINE__);
    e.SetLocation(ITK_LOCATION);
    e.SetDescription("Requested region is (at least partially) outside the largest possible region.");
    e.SetDataObject(inputPtr);
    throw e;
    }
}


template< class TInputImage, class TOutputImage, class TMaskPixelType>
const typename HaralickTextureFeaturesImageFilter< TInputImage, TOutputImage, TMaskPixelType>::InputImagePointer
HaralickTextureFeaturesImageFilter< TInputImage, TOutputImage, TMaskPixelType> 
::GetEnergyImage()
{
   typename TOutputImage::Pointer outputImage = const_cast< TOutputImage*> (this->GetOutput());
  /* typedef VectorIndexSelectionCastImageFilter< TOutputImage, TInputImage> VectorToScalarFilterType;
   typename VectorToScalarFilterType::Pointer indexFilter = VectorToScalarFilterType::New();
   indexFilter->SetInput(outputImage);
   indexFilter->SetIndex(0);
   indexFilter->Update();
   
   typename TInputImage::Pointer energyImage = TInputImage::New();
   energyImage->Graft(indexFilter->GetOutput());
  */
   typename TInputImage::Pointer energyImage = TInputImage::New();
   energyImage->CopyInformation(outputImage);
   energyImage->SetBufferedRegion(outputImage->GetBufferedRegion());
   energyImage->Allocate();
   energyImage->FillBuffer(0);
   ImageRegionIterator< TOutputImage > out(outputImage, outputImage->GetBufferedRegion());
   ImageRegionIterator< TInputImage > en(energyImage, energyImage->GetBufferedRegion());
   for(out.GoToBegin(), en.GoToBegin(); !out.IsAtEnd(); ++out, ++en)
   {
      typename TOutputImage::PixelType vec = out.Get();
      en.Set(vec[0]);
   }
   return energyImage;
}


template< class TInputImage, class TOutputImage, class TMaskPixelType>
const typename HaralickTextureFeaturesImageFilter< TInputImage, TOutputImage, TMaskPixelType>::InputImagePointer
HaralickTextureFeaturesImageFilter< TInputImage, TOutputImage, TMaskPixelType> 
::GetEntropyImage()
{
   typename TOutputImage::Pointer outputImage = const_cast< TOutputImage*> (this->GetOutput());
   /*typedef VectorIndexSelectionCastImageFilter< TOutputImage, TInputImage> VectorToScalarFilterType;
   typename VectorToScalarFilterType::Pointer indexFilter = VectorToScalarFilterType::New();
   indexFilter->SetInput(outputImage);
   indexFilter->SetIndex(1);
   indexFilter->Update();
   
   typename TInputImage::Pointer entropyImage = TInputImage::New();
   entropyImage->Graft(indexFilter->GetOutput());
   return entropyImage;
   */
   typename TInputImage::Pointer entropyImage = TInputImage::New();
   entropyImage->CopyInformation(outputImage);
   entropyImage->SetBufferedRegion(outputImage->GetBufferedRegion());
   entropyImage->Allocate();
   entropyImage->FillBuffer(0);
   ImageRegionIterator< TOutputImage > out(outputImage, outputImage->GetBufferedRegion());
   ImageRegionIterator< TInputImage > en(entropyImage, entropyImage->GetBufferedRegion());
   for(out.GoToBegin(), en.GoToBegin(); !out.IsAtEnd(); ++out, ++en)
   {
      typename TOutputImage::PixelType vec = out.Get();
      en.Set(vec[1]);
   }
   return entropyImage;
}


template< class TInputImage, class TOutputImage, class TMaskPixelType>
const typename HaralickTextureFeaturesImageFilter< TInputImage, TOutputImage, TMaskPixelType>::InputImagePointer
HaralickTextureFeaturesImageFilter< TInputImage, TOutputImage, TMaskPixelType> 
::GetCorrelationImage()
{
   typename TOutputImage::Pointer outputImage = const_cast< TOutputImage*> (this->GetOutput());
  /*   typedef VectorIndexSelectionCastImageFilter< TOutputImage, TInputImage> VectorToScalarFilterType;
   typename VectorToScalarFilterType::Pointer indexFilter = VectorToScalarFilterType::New();
   indexFilter->SetInput(outputImage);
   indexFilter->SetIndex(2);
   indexFilter->Update();
   
   typename TInputImage::Pointer correlationImage = TInputImage::New();
   correlationImage->Graft(indexFilter->GetOutput());
   return correlationImage;
  */
  typename TInputImage::Pointer correlationImage = TInputImage::New();
  correlationImage->CopyInformation(outputImage);
  correlationImage->SetBufferedRegion(outputImage->GetBufferedRegion());
  correlationImage->Allocate();
  correlationImage->FillBuffer(0);
  ImageRegionIterator< TOutputImage > out(outputImage, outputImage->GetBufferedRegion());
  ImageRegionIterator< TInputImage > en(correlationImage, correlationImage->GetBufferedRegion());
  for(out.GoToBegin(), en.GoToBegin(); !out.IsAtEnd(); ++out, ++en)
  {
     typename TOutputImage::PixelType vec = out.Get();
     en.Set(vec[2]);
  }
  return correlationImage;
}


template< class TInputImage, class TOutputImage, class TMaskPixelType>
const typename HaralickTextureFeaturesImageFilter< TInputImage, TOutputImage, TMaskPixelType>::InputImagePointer
HaralickTextureFeaturesImageFilter< TInputImage, TOutputImage, TMaskPixelType> 
::GetDifferenceMomentImage()
{
   typename TOutputImage::Pointer outputImage = const_cast< TOutputImage*> (this->GetOutput());
   /*typedef VectorIndexSelectionCastImageFilter< TOutputImage, TInputImage> VectorToScalarFilterType;
   typename VectorToScalarFilterType::Pointer indexFilter = VectorToScalarFilterType::New();
   indexFilter->SetInput(outputImage);
   indexFilter->SetIndex(3);
   indexFilter->Update();
   
   typename TInputImage::Pointer differenceMomentImage = TInputImage::New();
   differenceMomentImage->Graft(indexFilter->GetOutput());
   return differenceMomentImage;
   */
   typename TInputImage::Pointer differenceMomentImage = TInputImage::New();
   differenceMomentImage->CopyInformation(outputImage);
   differenceMomentImage->SetBufferedRegion(outputImage->GetBufferedRegion());
   differenceMomentImage->Allocate();
   differenceMomentImage->FillBuffer(0);
   ImageRegionIterator< TOutputImage > out(outputImage, outputImage->GetBufferedRegion());
   ImageRegionIterator< TInputImage > en(differenceMomentImage, differenceMomentImage->GetBufferedRegion());
   for(out.GoToBegin(), en.GoToBegin(); !out.IsAtEnd(); ++out, ++en)
   {
      typename TOutputImage::PixelType vec = out.Get();
      en.Set(vec[3]);
   }
   return differenceMomentImage;
}



template< class TInputImage, class TOutputImage, class TMaskPixelType>
const typename HaralickTextureFeaturesImageFilter< TInputImage, TOutputImage, TMaskPixelType>::
InputImagePointer HaralickTextureFeaturesImageFilter< TInputImage, TOutputImage, TMaskPixelType>
::GetInertiaImage()
{
   typename TOutputImage::Pointer outputImage = const_cast< TOutputImage*> (this->GetOutput());
   /*typedef VectorIndexSelectionCastImageFilter< TOutputImage, TInputImage> VectorToScalarFilterType;
   typename VectorToScalarFilterType::Pointer indexFilter = VectorToScalarFilterType::New();
   indexFilter->SetInput(outputImage);
   indexFilter->SetIndex(4);
   indexFilter->Update();
   
   typename TInputImage::Pointer inertiaImage = TInputImage::New();
   inertiaImage->Graft(indexFilter->GetOutput());
   return inertiaImage;
   */
   typename TInputImage::Pointer inertiaImage = TInputImage::New();
   inertiaImage->CopyInformation(outputImage);
   inertiaImage->SetBufferedRegion(outputImage->GetBufferedRegion());
   inertiaImage->Allocate();
   inertiaImage->FillBuffer(0);
   ImageRegionIterator< TOutputImage > out(outputImage, outputImage->GetBufferedRegion());
   ImageRegionIterator< TInputImage > en(inertiaImage, inertiaImage->GetBufferedRegion());
   for(out.GoToBegin(), en.GoToBegin(); !out.IsAtEnd(); ++out, ++en)
   {
      typename TOutputImage::PixelType vec = out.Get();
      en.Set(vec[4]);
   }
   return inertiaImage;
}


template< class TInputImage, class TOutputImage, class TMaskPixelType>
const typename HaralickTextureFeaturesImageFilter< TInputImage, TOutputImage, TMaskPixelType>::InputImagePointer
HaralickTextureFeaturesImageFilter< TInputImage, TOutputImage, TMaskPixelType> 
::GetClusterShadeImage()
{
   typename TOutputImage::Pointer outputImage = const_cast< TOutputImage*> (this->GetOutput());
   /*typedef VectorIndexSelectionCastImageFilter< TOutputImage, TInputImage> VectorToScalarFilterType;
   typename VectorToScalarFilterType::Pointer indexFilter = VectorToScalarFilterType::New();
   indexFilter->SetInput(outputImage);
   indexFilter->SetIndex(5);
   indexFilter->Update();
   
   typename TInputImage::Pointer clusterShadeImage = TInputImage::New();
   clusterShadeImage->Graft(indexFilter->GetOutput());
   return clusterShadeImage;
  */
  typename TInputImage::Pointer clusterShadeImage = TInputImage::New();
  clusterShadeImage->CopyInformation(outputImage);
  clusterShadeImage->SetBufferedRegion(outputImage->GetBufferedRegion());
  clusterShadeImage->Allocate();
  clusterShadeImage->FillBuffer(0);
  ImageRegionIterator< TOutputImage > out(outputImage, outputImage->GetBufferedRegion());
  ImageRegionIterator< TInputImage > en(clusterShadeImage, clusterShadeImage->GetBufferedRegion());
  for(out.GoToBegin(), en.GoToBegin(); !out.IsAtEnd(); ++out, ++en)
  {
     typename TOutputImage::PixelType vec = out.Get();
     en.Set(vec[5]);
  }
  return clusterShadeImage;
}


template< class TInputImage, class TOutputImage, class TMaskPixelType>
const typename HaralickTextureFeaturesImageFilter< TInputImage, TOutputImage, TMaskPixelType>::InputImagePointer
HaralickTextureFeaturesImageFilter< TInputImage, TOutputImage, TMaskPixelType> 
::GetClusterProminenceImage()
{
   typename TOutputImage::Pointer outputImage = const_cast< TOutputImage*> (this->GetOutput());
   /*typedef VectorIndexSelectionCastImageFilter< TOutputImage, TInputImage> VectorToScalarFilterType;
   typename VectorToScalarFilterType::Pointer indexFilter = VectorToScalarFilterType::New();
   indexFilter->SetInput(outputImage);
   indexFilter->SetIndex(6);
   indexFilter->Update();
   
   typename TInputImage::Pointer clusterProminenceImage = TInputImage::New();
   clusterProminenceImage->Graft(indexFilter->GetOutput());
   return clusterProminenceImage;
  */
  typename TInputImage::Pointer clusterPromImage = TInputImage::New();
  clusterPromImage->CopyInformation(outputImage);
  clusterPromImage->SetBufferedRegion(outputImage->GetBufferedRegion());
  clusterPromImage->Allocate();
  clusterPromImage->FillBuffer(0);
  ImageRegionIterator< TOutputImage > out(outputImage, outputImage->GetBufferedRegion());
  ImageRegionIterator< TInputImage > en(clusterPromImage, clusterPromImage->GetBufferedRegion());
  for(out.GoToBegin(), en.GoToBegin(); !out.IsAtEnd(); ++out, ++en)
  {
     typename TOutputImage::PixelType vec = out.Get();
     en.Set(vec[6]);
  }
  return clusterPromImage;

}


template< class TInputImage, class TOutputImage, class TMaskPixelType> 
void 
HaralickTextureFeaturesImageFilter< TInputImage, TOutputImage, TMaskPixelType>
::BeforeThreadedGenerateData()
{
   // initialize the outputImage
  typename TOutputImage::Pointer output = this->GetOutput();
  if(output->GetBufferedRegion().GetNumberOfPixels() == 0)
  {
     output->SetBufferedRegion( this->GetInput()->GetBufferedRegion());
     output->Allocate();
     typename TOutputImage::PixelType pix;
     for(unsigned i = 0; i < pix.GetVectorDimension(); i++)
     {
        pix[i] = 0;
     }
     output->FillBuffer(pix);
  } 
  
  m_Radius.Fill(m_WindowSize);
}


template< class TInputImage, class TOutputImage, class TMaskPixelType>
void 
HaralickTextureFeaturesImageFilter< TInputImage, TOutputImage, TMaskPixelType>
::AfterThreadedGenerateData()
{

    typename TOutputImage::Pointer outputImage = this->GetOutput();//const_cast< TOutputImage*> (this->GetOutput());

    ImageRegionIterator< TOutputImage > out(outputImage, outputImage->GetBufferedRegion());

    // write out the texture images 
    typedef Image< double, 3 > InternalImageType;
    typename InternalImageType::Pointer inImage = InternalImageType::New();
    inImage->CopyInformation(outputImage);
    inImage->SetBufferedRegion(outputImage->GetBufferedRegion());
    inImage->Allocate();
    inImage->FillBuffer(0);
   
   typedef RescaleIntensityImageFilter< InternalImageType, 
           InternalImageType > RescaleIntensityFilterType;

    for (unsigned n = 0; n < 7; ++n)
    {
       ImageRegionIterator< InternalImageType > in(inImage, inImage->GetBufferedRegion());
       for (in.GoToBegin(), out.GoToBegin(); !out.IsAtEnd(); ++in, ++out)
       {
           typename TOutputImage::PixelType vec = out.Get();
	   in.Set(vec[n]);
       }   

       if(m_NormalizeTextures)
       {
          std::cout<<" Normalizing features ... max "<<m_NormalizeMax<<" min "<<m_NormalizeMin<<std::endl;
          //typename TInputImage::Pointer scalarImage = TInputImage::New();
          //typedef VectorIndexSelectionCastImageFilter< TOutputImage, TInputImage> VectorToScalarFilterType;
          //for (unsigned n = 0; n < 7; ++n)
          //{
          //  typename VectorToScalarFilterType::Pointer vecToScalarFilter = VectorToScalarFilterType::New();
           // vecToScalarFilter->SetInput(outputImage);
           // vecToScalarFilter->SetIndex(n);
           // vecToScalarFilter->Update();

           // compute median filtering if that has been turned on
           if(m_ApplyMedianFiltering)
           {
              typedef MedianImageFilter< InternalImageType, InternalImageType > MedianImageFilterType;
	      typename MedianImageFilterType::Pointer medianFilter = MedianImageFilterType::New();
	      medianFilter->SetInput( inImage);
	  
	      typename TInputImage::SizeType radiusMedian;
	      radiusMedian.Fill(1);
	      medianFilter->SetRadius(radiusMedian);

	      medianFilter->Update();
	    
              // now rescale the intensity of the median filtered image
              typename RescaleIntensityFilterType::Pointer rescaleFilter = RescaleIntensityFilterType::New();
	      rescaleFilter->SetInput(medianFilter->GetOutput());
	      rescaleFilter->SetOutputMinimum( static_cast< unsigned long > (m_NormalizeMin) );
	      rescaleFilter->SetOutputMaximum( static_cast< unsigned long > (m_NormalizeMax ));
	      rescaleFilter->Update();

	      // put this back into the output
	      ImageRegionIterator< InternalImageType > scalar( rescaleFilter->GetOutput(), 
	            rescaleFilter->GetOutput()->GetBufferedRegion());
              ImageRegionIterator< TOutputImage > vector( outputImage, 
	            outputImage->GetBufferedRegion());
	      for(scalar.GoToBegin(), vector.GoToBegin(); !scalar.IsAtEnd(); 
	           ++scalar, ++vector)
	      {
                  typename TOutputImage::PixelType vec = vector.Get();
	          vec[n] = scalar.Get();
	          vector.Set(vec);
	      }
           }
           else
           {
             // simply rescale the intensity
              typename RescaleIntensityFilterType::Pointer rescaleFilter = RescaleIntensityFilterType::New();
	      rescaleFilter->SetInput(inImage);
	      rescaleFilter->SetOutputMinimum( m_NormalizeMin );
	      rescaleFilter->SetOutputMaximum( m_NormalizeMax );
	      rescaleFilter->Update();

	      // put this back into the output
	      ImageRegionIterator< InternalImageType > scalar( rescaleFilter->GetOutput(), 
	            rescaleFilter->GetOutput()->GetBufferedRegion());
              ImageRegionIterator< TOutputImage > vector( outputImage, 
	            outputImage->GetBufferedRegion());
	      for(scalar.GoToBegin(), vector.GoToBegin(); !scalar.IsAtEnd(); 
	           ++scalar, ++vector)
	      {
                  typename TOutputImage::PixelType vec = vector.Get();
	          vec[n] = scalar.Get();
	          vector.Set(vec);
	      }
           }
        //}
       }
        out.GoToBegin();
        for (in.GoToBegin(); !out.IsAtEnd(); ++in, ++out)
        { 
          typename TOutputImage::PixelType vec = out.Get();
	  in.Set(vec[n]);
        }
        typedef ImageFileWriter< InternalImageType > FileWriterType;
        typename FileWriterType::Pointer writer = FileWriterType::New();
        writer->SetInput(inImage);
        char str[128];
        sprintf(str, "HaralickFeature%d.mhd",n);
        writer->SetFileName(str);
        writer->Update();
    }
    
    // compute the texture statistics for the mask
    m_MeanEnergy = 0.0;
    m_MeanEntropy = 0.0;
    m_MeanCorrelation = 0.0;
    m_MeanDifferenceMoment = 0.0;
    m_MeanInertia = 0.0;
    m_MeanClusterShade = 0.0;
    m_MeanClusterProminence = 0.0;

    m_StdEnergy = 0.0;
    m_StdEntropy = 0.0;
    m_StdCorrelation = 0.0;
    m_StdDifferenceMoment = 0.0;
    m_StdInertia = 0.0;
    m_StdClusterShade = 0.0;
    m_StdClusterProminence = 0.0;

    m_KurtosisEnergy = 0.0;
    m_KurtosisEntropy = 0.0;
    m_KurtosisCorrelation = 0.0;
    m_KurtosisDifferenceMoment = 0.0;
    m_KurtosisInertia = 0.0;
    m_KurtosisClusterShade = 0.0;
    m_KurtosisClusterProminence = 0.0;

    ImageRegionIterator< MaskImageType > mask(m_Mask, m_Mask->GetBufferedRegion());
    //ImageRegionIterator< TOutputImage > out(outputImage, outputImage->GetBufferedRegion());

    std::vector< typename TOutputImage::PixelType > texVals;
    for (mask.GoToBegin(), out.GoToBegin(); !mask.IsAtEnd(); ++mask)
    {
      if(mask.Get())
      {
          typename TOutputImage::PixelType vec = out.Get();
	  texVals.push_back(vec);             
      }
    } 
    /*
    // compute the individual stats
    std::cout<<" size of texVals "<<texVals.size()<<std::endl;
    for (unsigned i = 0; i < texVals.size(); ++i)
    {
        m_MeanEnergy = m_MeanEnergy + texVals[i][0];
        m_MeanEntropy = m_MeanEntropy + texVals[i][1];
	m_MeanCorrelation = m_MeanCorrelation + texVals[i][2];
	m_MeanDifferenceMoment = m_MeanDifferenceMoment + texVals[i][3];
	m_MeanInertia = m_MeanInertia + texVals[i][4];
	m_MeanClusterShade = m_MeanClusterShade + texVals[i][5];
	m_MeanClusterProminence = m_MeanClusterProminence + texVals[i][6];
    }
    m_MeanEnergy /= static_cast< float > (texVals.size());
    m_MeanEntropy /= static_cast< float > (texVals.size());
    m_MeanCorrelation /= static_cast< float > (texVals.size());
    m_MeanDifferenceMoment /= static_cast< float > (texVals.size());
    m_MeanInertia /= static_cast< float > (texVals.size());
    m_MeanClusterShade /= static_cast< float > (texVals.size());
    m_MeanClusterProminence /= static_cast< float > (texVals.size());
    */
    //std::cout<<" Means Energy "<<m_MeanEnergy<<" Entropy "<<m_MeanEntropy<<" Correlation "<<m_MeanCorrelation<<" diffMoment "<<m_MeanDifferenceMoment<<" inertia "<<m_MeanInertia<<" shade "<<m_MeanClusterShade<<" prominence "<<m_MeanClusterProminence<<std::endl;
    

}


template< class TInputImage, class TOutputImage, class TMaskPixelType>
void
HaralickTextureFeaturesImageFilter< TInputImage, TOutputImage, TMaskPixelType>
::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
                       int threadId)
//::GenerateData()
{

//  Initialize();
  typename TInputImage::Pointer inputImage = const_cast< TInputImage * >(this->GetInput());
  
  InputPixelType minVal = NumericTraits< InputPixelType> ::max(InputPixelType() );
  InputPixelType maxVal = NumericTraits< InputPixelType > ::min(InputPixelType());
  ImageRegionIterator< TInputImage > it(  inputImage, outputRegionForThread); //inputImage->GetBufferedRegion()); //outputRegionForThread);
  for (it.GoToBegin(); !it.IsAtEnd(); ++it)
  {
     if(it.Get() < minVal)
     {
        minVal = it.Get();
     }
     if(it.Get() > maxVal)
     {
       maxVal = it.Get();
     }
  }
  //std::cout<<" Computing Haralick Texture "<<minVal<<" "<<maxVal<<" nbins "<<m_NumberOfBins<<std::endl;
  ConstantBoundaryCondition<MaskImageType> nbc;
  nbc.SetConstant(0);
  
  ConstNeighborhoodIterator<MaskImageType> bit;
  NeighborhoodIterator< TOutputImage > mit;
  
  //ImageRegionIterator<OutputImageType> it;
  
  typename OutputImageType::Pointer output = this->GetOutput();
  
  // Find the data-set boundary "faces"
  typename NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InputImageType>::FaceListType faceList;
  NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InputImageType> bC;
  faceList = bC(inputImage, outputRegionForThread, m_Radius); //inputImage->GetBufferedRegion(), m_Radius); //outputRegionForThread, m_Radius);

  typename NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InputImageType>::FaceListType::iterator fit;
  
  typedef Statistics::ScalarImageToTextureFeaturesFilter< TInputImage > TextureFilterType;
  
  // support progress methods/callbacks
  //ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels());
 
  typename TInputImage::SizeType imsize = inputImage->GetBufferedRegion().GetSize();
  //typename TInputImage::SizeType imsize = outputRegionForThread.GetSize();

  typename InputImageType::RegionType roiRegion;
  typename InputImageType::IndexType roiStart;
  typename InputImageType::SizeType roiSize;
  typename InputImageType::SpacingType roiSpacing = inputImage->GetSpacing();
  typename InputImageType::PointType roiOrigin = inputImage->GetOrigin();

  // Process each of the boundary faces.  These are N-d regions which border
  // the edge of the buffer.
  typedef ExtractImageFilter< InputImageType, InputImageType > ExtractFilterType;
  //typename ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
  //extractFilter->SetInput(inputImage);

  for (fit=faceList.begin(); fit != faceList.end(); ++fit)
    { 
    if(m_Mask.IsNotNull())
      {
        bit = ConstNeighborhoodIterator<MaskImageType>(m_Radius,
                                                        m_Mask, *fit);
        unsigned int neighborhoodSize = bit.Size();
        bit.OverrideBoundaryCondition(&nbc);
        bit.GoToBegin();                
        
        mit = NeighborhoodIterator< OutputImageType >(m_Radius, output, *fit); 
        mit.GoToBegin();

         
        //typedef RegionOfInterestImageFilter< InputImageType, InputImageType > ROIFilterType;
        //typename ROIFilterType::Pointer roiFilter = ROIFilterType::New();
     

        //typename InputImageType::RegionType iRegion;
        //typename InputImageType::IndexType iStart;
        //typename InputImageType::SizeType iSize;
        
        while(!bit.IsAtEnd())
        {
           if(bit.GetCenterPixel() != 0)
           {

              typename ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
              extractFilter->SetInput(inputImage);
              
              typename MaskImageType::IndexType first = bit.GetIndex(0);
              typename MaskImageType::IndexType last = bit.GetIndex(neighborhoodSize-1);
              for(unsigned n = 0; n < inputImage->GetImageDimension(); n++)
              {
                 roiStart[n] = first[n] >= 0 ? first[n] : 0;
                 roiSize[n] = (last[n] < imsize[n]) ? (last[n]-roiStart[n]) : (imsize[n]-1 - roiStart[n]);
                 //roiSize[n] += 1;
              }
              roiRegion.SetSize( roiSize );
              roiRegion.SetIndex( roiStart );
          
              extractFilter->SetExtractionRegion(roiRegion);
              try
              {
                 extractFilter->Update();
              }
              catch(itk::ExceptionObject &eExtract)
              {
                std::cout<<" Error in extract object "<<roiStart<<" "<<roiSize<<" "<<eExtract<<std::endl;
                typename TOutputImage::PixelType vPix;
                for (unsigned j = 0; j < vPix.GetVectorDimension(); j++)
                {
                   vPix[j] = 0;
                }
                mit.SetCenterPixel(vPix);
                ++bit;
                ++mit;
                //progress.CompletedPixel();
                continue;
              }    
             // roiFilter->SetRegionOfInterest( iRegion );
             // roiFilter->SetInput( inputImage );
            //  try
            //  {
            //     roiFilter->Update();
            //  }
            //  catch(itk::ExceptionObject &eroi)
            //  {
            //     std::cout<<" Error in roi Filter "<<iStart<<" "<<iSize<<"  "<<eroi<<std::endl;
            //     return;
            //  }
              InputPixelType minVal = NumericTraits< InputPixelType> ::max(InputPixelType() );
              InputPixelType maxVal = NumericTraits< InputPixelType > ::min(InputPixelType());
              ImageRegionIterator< TInputImage > it(  extractFilter->GetOutput(), extractFilter->GetOutput()->GetBufferedRegion()); //outputRegionForThread);
              for (it.GoToBegin(); !it.IsAtEnd(); ++it)
              {
                if(it.Get() < minVal)
                {
                  minVal = it.Get();
                }
                if(it.Get() > maxVal)
                {
                  maxVal = it.Get();
                }
              }
              //std::cout<<" max "<<maxVal<<" "<<minVal<<std::endl;
              typename TextureFilterType::Pointer textureFilter = TextureFilterType::New();
            //  textureFilter->SetInput( roiFilter->GetOutput());
                textureFilter->SetInput( extractFilter->GetOutput() );
                textureFilter->SetPixelValueMinMax( minVal, maxVal );
                //textureFilter->FastCalculationsOn();
                textureFilter->SetNumberOfBinsPerAxis( m_NumberOfBins);
                try
                {
                  textureFilter->Update();
                }
                catch(itk::ExceptionObject &etext)
                {
                   std::cout<<" error in texture filter "<<etext<<std::endl;
                   return;
                }
                // now get the texture features and add them to the output
                typename TextureFilterType::FeatureValueVector::Pointer featureVector = TextureFilterType::FeatureValueVector::New();
                featureVector = textureFilter->GetFeatureMeans();
                typename TOutputImage::PixelType vPix;
                for (unsigned j = 0; j < vPix.GetVectorDimension(); j++)
                {
                   vPix[j] = 0;
                }
                unsigned nfeatures = vPix.GetVectorDimension() < featureVector->Size() ? vPix.GetVectorDimension() : 
                                     featureVector->Size();
                for(unsigned j = 0; j < nfeatures; j++)
                {
                  vPix[j] = featureVector->GetElement(j);
		  //std::cout<<" texture : "<<vPix[j]<<" ";
                }
		//std::cout<<""<<std::endl;
                mit.SetCenterPixel(vPix);
           }	
           ++bit;
           ++mit;
           //progress.CompletedPixel();
        } // end while
      }
      else
      {

        //bit = ConstNeighborhoodIterator<MaskImageType>(m_Radius,
        //                                                m_Mask, *fit);
        //bit.OverrideBoundaryCondition(&nbc);
        //bit.GoToBegin();

        mit = NeighborhoodIterator< OutputImageType >(m_Radius, output, *fit);
        mit.GoToBegin();

        unsigned int neighborhoodSize = mit.Size();
        typedef RegionOfInterestImageFilter< InputImageType, InputImageType > ROIFilterType;
        typename ROIFilterType::Pointer roiFilter = ROIFilterType::New();
        typename InputImageType::RegionType iRegion;
        typename InputImageType::IndexType iStart;
        typename InputImageType::SizeType iSize;


         while(!mit.IsAtEnd())
         {
           std::cout<<" MASK IS NOT SET "<<std::endl;
           typename MaskImageType::IndexType first = mit.GetIndex(0);
           typename MaskImageType::IndexType last = mit.GetIndex(neighborhoodSize-1);
           for(unsigned n = 0; n < inputImage->GetImageDimension(); n++)
           {
              iStart[n] = first[n] >= 0 ? first[n] : 0;
              iSize[n] = (last[n] < imsize[n]) ? (last[n]-iStart[n]) : (imsize[n]-1 - iStart[n]);
              iSize[n] += 1;
           }
           iRegion.SetSize( iSize );
           iRegion.SetIndex( iStart );

           roiFilter->SetRegionOfInterest( iRegion );
           roiFilter->SetInput( inputImage );
           roiFilter->Update();

           typename TextureFilterType::Pointer textureFilter = TextureFilterType::New();
           textureFilter->SetInput( roiFilter->GetOutput());
           textureFilter->SetPixelValueMinMax( minVal, maxVal );
           //textureFilter->FastCalculationsOn();
           textureFilter->SetNumberOfBinsPerAxis( m_NumberOfBins);
           textureFilter->Update();

           // now get the texture features and add them to the output
           typename TextureFilterType::FeatureValueVector::Pointer featureVector = TextureFilterType::FeatureValueVector::New();
           featureVector = textureFilter->GetFeatureMeans();
           typename TOutputImage::PixelType vPix;
           for (unsigned j = 0; j < vPix.Size(); j++)
           {
              vPix[j] = 0;
           }
           unsigned nfeatures = vPix.Size() < featureVector->Size() ? vPix.Size() :
                                featureVector->Size();
           for(unsigned j = 0; j < nfeatures; j++)
           {
              vPix[j] = featureVector->GetElement(j);
           }
           mit.SetCenterPixel(vPix);
           //++bit;
           ++mit;
           //progress.CompletedPixel();
         }        
      }
    }
}

/**
 * Standard "PrintSelf" method
 */
template <class TInputImage, class TOutput, class TMaskPixelType>
void
HaralickTextureFeaturesImageFilter<TInputImage, TOutput, TMaskPixelType>
::PrintSelf(
  std::ostream& os, 
  Indent indent) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "Radius: " << m_Radius << std::endl;
  os << indent << "Number of Bins: " << m_NumberOfBins << std::endl;
}

} // end namespace itk

#endif
