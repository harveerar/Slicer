/**********************************************************************
 * vtkItkHaralickTextureFeaturesImageFilter
 * Implements wrapper for the itkHaralickTextureFeaturesImageFilter
 * This implemnents n-class segmentation
 **********************************************************************/

#include "vtkITKHaralickTextureFeaturesImageFilter.h"

// VTK includes
#include <vtkImageCast.h>
#include <vtkImageData.h>
#include <vtkObjectFactory.h>

// ITK includes
#include <itkHaralickTextureFeaturesImageFilter.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageRegionIterator.h>
#include <itkRegionOfInterestImageFilter.h>
#include <itkVectorImage.h>
#include <itkVector.h>

#include <itkListSample.h>
#include <itkSampleToHistogramFilter.h>
#include <itkHistogram.h>

#include <vector>
#include <algorithm>

//-----------------------------------------------------------------------------
vtkCxxRevisionMacro(vtkITKHaralickTextureFeaturesImageFilter, "$Revision: 1.1 $");
vtkStandardNewMacro(vtkITKHaralickTextureFeaturesImageFilter);

//-----------------------------------------------------------------------------
////////// These types are not defined in itk ////////////
#ifdef vtkTemplateMacroCase_ui64
#undef vtkTemplateMacroCase_ui64
# define vtkTemplateMacroCase_ui64(typeN, type, call)
#endif
#ifdef vtkTemplateMacroCase_si64
#undef vtkTemplateMacroCase_si64
# define vtkTemplateMacroCase_si64(typeN, type, call)
#endif
#ifdef vtkTemplateMacroCase_ll
#undef vtkTemplateMacroCase_ll
# define vtkTemplateMacroCase_ll(typeN, type, call)
#endif

//-----------------------------------------------------------------------------
// Local Function: not method.
void vtkITKImageHaralickTextureFeatureHandleProgressEvent(itk::Object *caller,
                                           const itk::EventObject& vtkNotUsed(eventObject),
                                           void *clientdata)
{

  itk::ProcessObject *itkFilter = static_cast<itk::ProcessObject*>(caller);
  vtkProcessObject *vtkFilter = static_cast<vtkProcessObject*>(clientdata);
  if (itkFilter && vtkFilter )
    {
    vtkFilter->UpdateProgress( itkFilter->GetProgress() );
    }
};


//-----------------------------------------------------------------------------
//// 3D filter
template<class IT1, class OT>
void vtkITKImageHaralickTextureFeaturesExecute3D(vtkImageData *inData,
  IT1 *inPtr1, OT *inPtr2,
  IT1 *outputEnergy, IT1* outputEntropy,
  IT1 *outputCorrelation, IT1* outputDiffMoment,
  IT1 *outputInertia, IT1 *outputClusterShade, 
  IT1 *outputClusterProminence, 
  double &Radius,
  double &NumberOfBins,
  bool NormalizeTextures,
  bool ApplyMedianFilter,
  float &NormMaximum,
  float &NormMinimum, 
  std::vector< int > &labelIndexMap,
  std::vector< float > &meanEnergy, std::vector< float > &meanEntropy,
  std::vector< float > &meanCorrelation, std::vector< float > &meanInertia,
  std::vector< float > &meanDifferenceMoment, std::vector< float > &meanClusterShade,
  std::vector< float > &meanClusterProminence, std::vector< float > &stdEnergy,
  std::vector< float > &stdEntropy, std::vector< float > &stdCorrelation,
  std::vector< float > &stdDifferenceMoment, std::vector< float > &stdInertia,
  std::vector< float > &stdClusterShade, std::vector< float > &stdClusterProminence,
  std::vector< float > &kurtosisEnergy, std::vector< float > &kurtosisEntropy,
  std::vector< float > &kurtosisCorrelation, std::vector< float > &kurtosisDifferenceMoment,
  std::vector< float > &kurtosisInertia, std::vector< float > &kurtosisClusterShade,
  std::vector< float > &kurtosisClusterProminence,
  std::vector< float > &skewEnergy, std::vector< float > &skewEntropy,
  std::vector< float > &skewCorrelation, std::vector< float > &skewDifferenceMoment,
  std::vector< float > &skewInertia, std::vector< float > &skewClusterShade,
  std::vector< float > &skewClusterProminence, 
  std::vector< float > &medianEnergy, std::vector< float > &medianEntropy,
  std::vector< float > &medianCorrelation, std::vector< float > &medianDifferenceMoment,
  std::vector< float > &medianInertia, std::vector< float > &medianClusterShade, 
  std::vector< float > &medianClusterProminence
 )
{
  typedef itk::Image<IT1, 3> FeatureImageType;
  typename FeatureImageType::Pointer image = FeatureImageType::New();
  
  typedef itk::Image<OT, 3> LabelImageType;

  typename LabelImageType::Pointer labelImage = LabelImageType::New();

  typename LabelImageType::Pointer maskImage = LabelImageType::New();

  typedef itk::Vector< double, 7 > TexturePixelType;
  typedef itk::Image< TexturePixelType, 3> TextureImageType;

  typename TextureImageType::Pointer outputImage = TextureImageType::New();

  int dims[3];
  int extent[6];
  double spacing[3], origin[3];

  inData->GetDimensions(dims);
  inData->GetExtent(extent);
  inData->GetOrigin(origin);
  inData->GetSpacing(spacing);

  image->SetOrigin( origin );
  image->SetSpacing( spacing );

  typename FeatureImageType::RegionType region;
  typename FeatureImageType::IndexType index;
  typename FeatureImageType::SizeType size;
  index[0] = extent[0];
  index[1] = extent[2];
  index[2] = extent[4];
  region.SetIndex( index );
  size[0] = extent[1] - extent[0] + 1;
  size[1] = extent[3] - extent[2] + 1;
  size[2] = extent[5] - extent[4] + 1;
  region.SetSize( size );
  image->SetRegions(region);

  image->GetPixelContainer()->SetImportPointer(inPtr1, dims[0]*dims[1]*dims[2], false);

  labelImage->SetOrigin( origin );
  labelImage->SetSpacing( spacing );
  labelImage->SetRegions( region );
  labelImage->GetPixelContainer()->SetImportPointer(inPtr2, dims[0]*dims[1]*dims[2], false);

  maskImage->CopyInformation(image);
  maskImage->SetBufferedRegion( image->GetBufferedRegion() );
  maskImage->Allocate();
  maskImage->FillBuffer( 0 );

  itk::ImageRegionIteratorWithIndex< LabelImageType > label(labelImage, labelImage->GetBufferedRegion() );
  // get all the unique labels. Construct a separate ROI for each label and compute the Haralick Textures inside those 
  // regions
  std::vector< int > uniqueLabels;
  for (label.GoToBegin(); !label.IsAtEnd(); ++label)
  {
    if(label.Get() > 0)
    {
      bool found = false;
      for(unsigned n = 0; uniqueLabels.size() > 0 && !found && n < uniqueLabels.size(); ++n)
      {
	 found = uniqueLabels[n] == label.Get();
      }
      if(!found)
      {
	 uniqueLabels.push_back(label.Get());
      }
     }   
  }
  
  if(uniqueLabels.size() < 1)
  {
      std::cout<<" number of unique labels "<<uniqueLabels.size()<<std::endl;
      return;
  }

  labelIndexMap.resize(uniqueLabels.size());
  for (unsigned n = 0; n < labelIndexMap.size(); ++n)
  {
    labelIndexMap[n] = uniqueLabels[n];
  }
  std::cout<<" number of unique labels "<<labelIndexMap.size()<<std::endl;

  // Allocate the Texture feature images
  typename LabelImageType::Pointer energyImage = LabelImageType::New();
  typename LabelImageType::Pointer entropyImage = LabelImageType::New();
  typename LabelImageType::Pointer correlationImage = LabelImageType::New();
  typename LabelImageType::Pointer inertiaImage = LabelImageType::New();
  typename LabelImageType::Pointer differenceMomentImage = LabelImageType::New();
  typename LabelImageType::Pointer clusterShadeImage = LabelImageType::New();
  typename LabelImageType::Pointer clusterProminenceImage = LabelImageType::New();
 
  energyImage->CopyInformation(image);
  energyImage->SetBufferedRegion(image->GetBufferedRegion());
  energyImage->Allocate();
  energyImage->FillBuffer(0);

  entropyImage->CopyInformation(image);
  entropyImage->SetBufferedRegion(image->GetBufferedRegion());
  entropyImage->Allocate();
  entropyImage->FillBuffer(0);

  correlationImage->CopyInformation(image);
  correlationImage->SetBufferedRegion(image->GetBufferedRegion());
  correlationImage->Allocate();
  correlationImage->FillBuffer(0);

  inertiaImage->CopyInformation(image);
  inertiaImage->SetBufferedRegion(image->GetBufferedRegion());
  inertiaImage->Allocate();
  inertiaImage->FillBuffer(0);

  differenceMomentImage->CopyInformation(image);
  differenceMomentImage->SetBufferedRegion(image->GetBufferedRegion());
  differenceMomentImage->Allocate();
  differenceMomentImage->FillBuffer(0);

  clusterShadeImage->CopyInformation(image);
  clusterShadeImage->SetBufferedRegion(image->GetBufferedRegion());
  clusterShadeImage->Allocate();
  clusterShadeImage->FillBuffer(0);

  clusterProminenceImage->CopyInformation(image);
  clusterProminenceImage->SetBufferedRegion(image->GetBufferedRegion());
  clusterProminenceImage->Allocate();
  clusterProminenceImage->FillBuffer(0);

   unsigned int nLabels = uniqueLabels.size();
   meanEnergy.resize(nLabels);
   meanEntropy.resize(nLabels);
   meanCorrelation.resize(nLabels);
   meanDifferenceMoment.resize(nLabels);
   meanInertia.resize(nLabels);
   meanClusterShade.resize(nLabels);
   meanClusterProminence.resize(nLabels);

   stdEnergy.resize(nLabels);
   stdEntropy.resize(nLabels);
   stdCorrelation.resize(nLabels);
   stdDifferenceMoment.resize(nLabels);
   stdInertia.resize(nLabels);
   stdClusterShade.resize(nLabels);
   stdClusterProminence.resize(nLabels);
    
   kurtosisEnergy.resize(nLabels);
   kurtosisEntropy.resize(nLabels);
   kurtosisCorrelation.resize(nLabels);
   kurtosisDifferenceMoment.resize(nLabels);
   kurtosisInertia.resize(nLabels);
   kurtosisClusterShade.resize(nLabels);
   kurtosisClusterProminence.resize(nLabels); 

   skewEnergy.resize(nLabels);
   skewEntropy.resize(nLabels);
   skewCorrelation.resize(nLabels);
   skewDifferenceMoment.resize(nLabels);
   skewInertia.resize(nLabels);
   skewClusterShade.resize(nLabels);
   skewClusterProminence.resize(nLabels);

   medianEnergy.resize(nLabels);
   medianEntropy.resize(nLabels);
   medianCorrelation.resize(nLabels);
   medianDifferenceMoment.resize(nLabels);
   medianInertia.resize(nLabels);
   medianClusterShade.resize(nLabels);
   medianClusterProminence.resize(nLabels);

   for(unsigned nl = 0; nl < uniqueLabels.size(); ++nl)
   {
    // set up the ROI mask for processing using the Haralick Texture filter
    typename FeatureImageType::IndexType startROI;
    typename FeatureImageType::IndexType endROI;
    typename FeatureImageType::SizeType sizeROI;
    typename FeatureImageType::RegionType regionROI;

    startROI[0] = index[0];
    startROI[1] = index[1];
    startROI[2] = index[2];
  
    endROI[0] = 0;
    endROI[1] = 0;
    endROI[2] = 0;

    sizeROI[0] = size[0];
    sizeROI[1] = size[1];
    sizeROI[2] = size[2];
  
    bool setFirst = false;
    // Extract the ROI using the label mask
    for(label.GoToBegin(); !label.IsAtEnd(); ++label)
      {
	if(label.Get() == uniqueLabels[nl])
	  {
	    typename FeatureImageType::IndexType indx = label.GetIndex();
	    for (unsigned int n = 0; n < 3; ++n)
	      {
		if(!setFirst || startROI[n] >  indx[n])
		  {
		    startROI[n] = indx[n];
		  }
		if(!setFirst || endROI[n] < indx[n])
		  {
		    endROI[n] = indx[n];
		  }
	      }
	    if(!setFirst)
	      {
		setFirst = true;
	      }
	  }
      }

    for (unsigned n = 0; n < 3; ++n)
      {
	sizeROI[n] = endROI[n] - startROI[n]; 
	if(sizeROI[n]+Radius/2 < size[n])
	  {
	    sizeROI[n] += Radius/2;
	  }
	if(startROI[n]-Radius/2 >= index[n])
	  {
	    startROI[n] -= Radius/2;
	  }
      }

    std::cout<<" roi start "<<startROI[0]<<", "<<startROI[1]<<", "<<startROI[2]<<std::endl;
    std::cout<<" roi size "<<sizeROI[0]<<", "<<sizeROI[1]<<", "<<sizeROI[2]<<std::endl;

    regionROI.SetIndex(startROI);
    regionROI.SetSize(sizeROI);

    // fill mask with all ones
    itk::ImageRegionIterator< LabelImageType > mask(maskImage, regionROI );
    itk::ImageRegionIterator< LabelImageType > labelROI(labelImage, regionROI);
    for (mask.GoToBegin(), labelROI.GoToBegin(); !mask.IsAtEnd(); 
	 ++mask, ++labelROI)
      {
	if(labelROI.Get() == uniqueLabels[nl])
	{
	  mask.Set(255);
	}
	else
	{
	    mask.Set(0);
	}
      }

    // set up the Haralick Texture filter  
    typedef itk::HaralickTextureFeaturesImageFilter<FeatureImageType, TextureImageType, OT> FilterType;
    typename FilterType::Pointer filter = FilterType::New();
    
    // filter->AddObserver(itk::ProgressEvent(), progressCommand );
    //typename FeatureImageType::SizeType radius;
    //radius.Fill(Radius);
    std::cout<<" Normalize textures value "<<NormalizeTextures<<std::endl;
    
    filter->SetInput(image);
    filter->SetMask(maskImage);
    filter->SetWindowSize(Radius);
    filter->SetNumberOfBins(NumberOfBins);
    filter->SetNormalizeTextures(NormalizeTextures);
    filter->SetApplyMedianFiltering(ApplyMedianFilter);
    filter->SetNormalizeMax(NormMaximum);
    filter->SetNormalizeMin(NormMinimum);

    std::cout<<" Now Running the TEXTURE FILTER ... "<<std::endl;
    filter->Update();
    std::cout<<" DONE running the texture analysis ... "<<std::endl;

    typename TextureImageType::Pointer texImage = TextureImageType::New();
    texImage->Graft(filter->GetOutput());
    itk::ImageRegionIterator< TextureImageType > tex(texImage, regionROI);
    itk::ImageRegionIterator< LabelImageType > energy( energyImage, regionROI);
    itk::ImageRegionIterator< LabelImageType > entropy( entropyImage, regionROI);
    itk::ImageRegionIterator< LabelImageType > correlation( correlationImage, regionROI);
    itk::ImageRegionIterator< LabelImageType > inertia( inertiaImage, regionROI);
    itk::ImageRegionIterator< LabelImageType > diffMoment( differenceMomentImage, 
							regionROI);
    itk::ImageRegionIterator< LabelImageType > clusterShade( clusterShadeImage, 
							  regionROI);
    itk::ImageRegionIterator< LabelImageType > clusterProminence( clusterProminenceImage, 
							       regionROI);
    tex.GoToBegin(); energy.GoToBegin();
    entropy.GoToBegin(); correlation.GoToBegin();
    diffMoment.GoToBegin(); inertia.GoToBegin(); 
    clusterShade.GoToBegin();
    clusterProminence.GoToBegin();
    mask.GoToBegin();

    std::vector< float > valsCorrelation;
    std::vector< float > valsEnergy;
    std::vector< float > valsEntropy;
    std::vector< float > valsDifferenceMoment;
    std::vector< float > valsInertia;
    std::vector< float > valsClusterShade;
    std::vector< float > valsClusterProminence;
    
    for (; !mask.IsAtEnd(); ++tex, ++energy,++mask, 
           ++entropy, ++correlation, ++diffMoment, ++inertia, 
           ++clusterShade, ++clusterProminence)
    {
	if(mask.Get())
	  {

	    TexturePixelType pix = tex.Get();
	    energy.Set(static_cast< OT >(pix[0]));
	    entropy.Set(static_cast< OT > (pix[1]));
	    correlation.Set(static_cast< OT >(pix[2]));
	    diffMoment.Set(static_cast< OT > (pix[3]));
	    inertia.Set(static_cast< OT > (pix[4]));
	    clusterShade.Set(static_cast< OT > (pix[5]));
	    clusterProminence.Set(static_cast< OT > (pix[6]));

	    valsEnergy.push_back(pix[0]);
	    valsEntropy.push_back(pix[1]);
	    valsCorrelation.push_back(pix[2]);
	    valsDifferenceMoment.push_back(pix[3]);
	    valsInertia.push_back(pix[4]);
	    valsClusterShade.push_back(pix[5]);
	    valsClusterProminence.push_back(pix[6]);

	    meanEnergy[nl] = meanEnergy[nl] + pix[0];
	    meanEntropy[nl] = meanEntropy[nl] + pix[1];
	    meanCorrelation[nl] = meanCorrelation[nl] + pix[2];
	    meanDifferenceMoment[nl] = meanDifferenceMoment[nl] + pix[3];
	    meanInertia[nl] = meanInertia[nl] + pix[4];
	    meanClusterShade[nl] = meanClusterShade[nl] + pix[5];
	    meanClusterProminence[nl] = meanClusterProminence[nl] + pix[6];
	  }
    }
 
    std::sort(valsEnergy.begin(), valsEnergy.end());
    std::sort(valsEntropy.begin(), valsEntropy.end());
    std::sort(valsCorrelation.begin(), valsCorrelation.end());
    std::sort(valsDifferenceMoment.begin(), valsDifferenceMoment.end());
    std::sort(valsInertia.begin(), valsInertia.end());
    std::sort(valsClusterShade.begin(), valsClusterShade.end());
    std::sort(valsClusterProminence.begin(), valsClusterProminence.end());

    // measure the mean, kurtosis and std
    meanEnergy[nl] = (valsEnergy.size() > 0) ? 
	meanEnergy[nl] / valsEnergy.size() : meanEnergy[nl];
    meanEntropy[nl] = (valsEntropy.size() > 0) ? 
	meanEntropy[nl] / valsEntropy.size() : meanEntropy[nl];
    meanCorrelation[nl] = (valsCorrelation.size() > 0) ? 
	meanCorrelation[nl] / valsCorrelation.size() : meanCorrelation[nl];
    meanDifferenceMoment[nl] = (valsDifferenceMoment.size() > 0) ? 
	meanDifferenceMoment[nl] / valsDifferenceMoment.size() : meanDifferenceMoment[nl];
    meanInertia[nl] = (valsInertia.size() > 0) ? 
	meanInertia[nl] / valsInertia.size() : meanInertia[nl];
    meanClusterShade[nl] = (valsClusterShade.size() > 0) ? 
	meanClusterShade[nl] / valsClusterShade.size() : meanClusterShade[nl];
    meanClusterProminence[nl] = (valsClusterProminence.size() > 0) ? 
	meanClusterProminence[nl] / valsClusterProminence.size() : meanClusterProminence[nl];

    int size = valsEnergy.size();
    int mid = size/2;

    stdEnergy[nl] = 0.0;
    kurtosisEnergy[nl] = 0.0;
    skewEnergy[nl] = 0.0;
    medianEnergy[nl] = valsEnergy[mid];

    stdEntropy[nl] = 0.0;
    kurtosisEntropy[nl] = 0.0;
    skewEntropy[nl] = 0.0;
    medianEntropy[nl] = valsEntropy[mid];

    stdCorrelation[nl] = 0.0;
    kurtosisCorrelation[nl] = 0.0;
    skewCorrelation[nl] = 0.0;
    medianCorrelation[nl] = valsCorrelation[mid];

    stdDifferenceMoment[nl] = 0.0;
    kurtosisDifferenceMoment[nl] = 0.0;
    skewDifferenceMoment[nl] = 0.0;
    medianDifferenceMoment[nl] = valsDifferenceMoment[mid];

    stdInertia[nl] = 0.0;
    kurtosisInertia[nl]= 0.0;
    skewInertia[nl] = 0.0;
    medianInertia[nl] = valsInertia[mid];

    stdClusterShade[nl] = 0.0;
    kurtosisClusterShade[nl] = 0.0;
    skewClusterShade[nl] = 0.0;
    medianClusterShade[nl] = valsClusterShade[mid];

    stdClusterProminence[nl] = 0.0;
    kurtosisClusterProminence[nl] = 0.0;
    skewClusterProminence[nl] = 0.0;
    medianClusterProminence[nl] = valsClusterProminence[mid];

    for (unsigned v = 0; v < valsEnergy.size(); ++v)
    {
       stdEnergy[nl] += (valsEnergy[v]-meanEnergy[nl])*(valsEnergy[v]-meanEnergy[nl]);
       kurtosisEnergy[nl] += std::pow((valsEnergy[v] - meanEnergy[nl]),4);
       skewEnergy[nl]+= std::pow((valsEnergy[v] - meanEnergy[nl]),3); 

       stdEntropy[nl] += (valsEntropy[v]-meanEntropy[nl])*(valsEntropy[v]-meanEntropy[nl]);
       kurtosisEntropy[nl] += std::pow((valsEntropy[v] - meanEntropy[nl]),4);
       skewEntropy[nl] += std::pow((valsEntropy[v] - meanEntropy[nl]),3);

       stdCorrelation[nl] += (valsCorrelation[v]-meanCorrelation[nl])*(valsCorrelation[v]-meanCorrelation[nl]);
       kurtosisCorrelation[nl] += std::pow((valsCorrelation[v] - meanCorrelation[nl]),4);
       skewCorrelation[nl] += std::pow((valsCorrelation[v] - meanCorrelation[nl]),3);

       stdDifferenceMoment[nl] += (valsDifferenceMoment[v]-meanDifferenceMoment[nl])*(valsDifferenceMoment[v]-meanDifferenceMoment[nl]);
       kurtosisDifferenceMoment[nl] += std::pow((valsDifferenceMoment[v] - meanDifferenceMoment[nl]),4);
       skewDifferenceMoment[nl] += std::pow((valsDifferenceMoment[v] - meanDifferenceMoment[nl]),3);


       stdInertia[nl] += (valsInertia[v]-meanInertia[nl])*(valsInertia[v]-meanInertia[nl]);
       kurtosisInertia[nl] += std::pow((valsInertia[v] - meanInertia[nl]),4);
       skewInertia[nl] += std::pow((valsInertia[v] - meanInertia[nl]),3);
	
       stdClusterShade[nl] += (valsClusterShade[v]-meanClusterShade[nl])*(valsClusterShade[v]-meanClusterShade[nl]);
       kurtosisClusterShade[nl] += std::pow((valsClusterShade[v] - meanClusterShade[nl]),4);
       skewClusterShade[nl] += std::pow((valsClusterShade[v] - meanClusterShade[nl]),3);

       stdClusterProminence[nl] += (valsClusterProminence[v]-meanClusterProminence[nl])*(valsClusterProminence[v]-meanClusterProminence[nl]);
       kurtosisClusterProminence[nl] += std::pow((valsClusterProminence[v] - meanClusterProminence[nl]),4);
       skewClusterProminence[nl] += std::pow((valsClusterProminence[v] - meanClusterProminence[nl]),3);
    }
    stdEnergy[nl] /= valsEnergy.size();
    stdEnergy[nl] = std::sqrt(stdEnergy[nl]);
    stdEnergy[nl] = stdEnergy[nl] > 0 ? stdEnergy[nl] : 1;
    kurtosisEnergy[nl] /= ((valsEnergy.size()-1)*std::pow(stdEnergy[nl],4));
    kurtosisEnergy[nl] -= 3;

    skewEnergy[nl] /= ((valsEnergy.size()-1)*std::pow(stdEnergy[nl],3));

    stdEntropy[nl] /= valsEntropy.size();
    stdEntropy[nl] = std::sqrt(stdEntropy[nl]);
    stdEntropy[nl] = stdEntropy[nl] > 0 ? stdEntropy[nl] : 1;
    kurtosisEntropy[nl] /= ((valsEntropy.size()-1)*std::pow(stdEntropy[nl],4));
    kurtosisEntropy[nl] -= 3;
    skewEntropy[nl] /= ((valsEntropy.size()-1)*std::pow(stdEntropy[nl],3));


    stdCorrelation[nl] /= valsCorrelation.size();
    stdCorrelation[nl] = std::sqrt(stdCorrelation[nl]);
    stdCorrelation[nl] = stdCorrelation[nl] > 0 ? stdCorrelation[nl] : 1;
    kurtosisCorrelation[nl] /= ((valsCorrelation.size()-1)*std::pow(stdCorrelation[nl],4));
    kurtosisCorrelation[nl] -= 3;
    skewCorrelation[nl] /= ((valsCorrelation.size()-1)*std::pow(stdCorrelation[nl],3));

    stdDifferenceMoment[nl] /= valsDifferenceMoment.size();
    stdDifferenceMoment[nl] = std::sqrt(stdDifferenceMoment[nl]);
    stdDifferenceMoment[nl] = stdDifferenceMoment[nl] > 0 ? stdDifferenceMoment[nl] : 1;
    kurtosisDifferenceMoment[nl] /= ((valsDifferenceMoment.size()-1)*std::pow(stdDifferenceMoment[nl],4));
    kurtosisDifferenceMoment[nl] -= 3;
    skewDifferenceMoment[nl] /= ((valsDifferenceMoment.size()-1)*std::pow(stdDifferenceMoment[nl],3));

    stdInertia[nl] /= valsInertia.size();
    stdInertia[nl] = std::sqrt(stdInertia[nl]);
    stdInertia[nl] = stdInertia[nl] > 0 ? stdInertia[nl] : 1;
    kurtosisInertia[nl] /= ((valsInertia.size()-1)*std::pow(stdInertia[nl],4));
    kurtosisInertia[nl] -= 3;
    skewInertia[nl] /= ((valsInertia.size()-1)*std::pow(stdInertia[nl],3));

    
    stdClusterShade[nl] = std::sqrt(stdClusterShade[nl]/valsClusterShade.size());
    stdClusterShade[nl] = stdClusterShade[nl] > 0 ? stdClusterShade[nl] : 1;
    kurtosisClusterShade[nl] /= ((valsClusterShade.size()-1)*std::pow(stdClusterShade[nl],4));
    kurtosisClusterShade[nl] -= 3;
    skewClusterShade[nl] /= ((valsClusterShade.size()-1)*std::pow(stdClusterShade[nl],3));

    stdClusterProminence[nl] = std::sqrt(stdClusterProminence[nl]/valsClusterProminence.size());
    stdClusterProminence[nl] = stdClusterProminence[nl] > 0 ? stdClusterProminence[nl] : 1;
    kurtosisClusterProminence[nl] /= ((valsClusterProminence.size()-1)*std::pow(stdClusterProminence[nl],4));
    kurtosisClusterProminence[nl] -= 3;
    skewClusterProminence[nl] /= ((valsClusterProminence.size()-1)*std::pow(stdClusterProminence[nl],3));

    std::cout<<" Label "<<nl<<" Energy mean "<<meanEnergy[nl]<<" kurtosis "<<kurtosisEnergy[nl]<<" std "<<stdEnergy[nl]<<std::endl;
    std::cout<<" Entropy mean "<<meanEntropy[nl]<<" kurtosis "<<kurtosisEntropy[nl]<<" std "<<stdEntropy[nl]<<std::endl;
    std::cout<<" Correlation mean "<<meanCorrelation[nl]<<" kurtosis "<<kurtosisCorrelation[nl]<<" std "<<stdCorrelation[nl]<<std::endl;
    std::cout<<" Diff Moment mean "<<meanDifferenceMoment[nl]<<" kurtosis "<<kurtosisDifferenceMoment[nl]<<" std "<<stdDifferenceMoment[nl]<<std::endl;
    std::cout<<" Inertia mean "<<meanInertia[nl]<<" kurtosis "<<kurtosisInertia[nl]<<" std "<<stdInertia[nl]<<std::endl;
    std::cout<<" ClusterShade mean "<<meanClusterShade[nl]<<" kurtosis "<<kurtosisClusterShade[nl]<<" std "<<stdClusterShade[nl]<<std::endl;
    std::cout<<" ClusterProminence mean "<<meanClusterProminence[nl]<<" kurtosis "<<kurtosisClusterProminence[nl]<<" std "<<stdClusterProminence[nl]<<std::endl;
  }

  // copy the individual textures
 
  memcpy(outputEnergy, energyImage->GetBufferPointer(),
         energyImage->GetBufferedRegion().GetNumberOfPixels()*sizeof(OT) );

  memcpy(outputEntropy, entropyImage->GetBufferPointer(),
         entropyImage->GetBufferedRegion().GetNumberOfPixels()*sizeof(OT) );

  memcpy(outputCorrelation, correlationImage->GetBufferPointer(),
         correlationImage->GetBufferedRegion().GetNumberOfPixels()*sizeof(OT) );

  memcpy(outputDiffMoment, differenceMomentImage->GetBufferPointer(),
         differenceMomentImage->GetBufferedRegion().GetNumberOfPixels()*sizeof(OT) );

  memcpy(outputInertia, inertiaImage->GetBufferPointer(),
         inertiaImage->GetBufferedRegion().GetNumberOfPixels()*sizeof(OT) );

  memcpy(outputClusterShade, clusterShadeImage->GetBufferPointer(),
         clusterShadeImage->GetBufferedRegion().GetNumberOfPixels()*sizeof(OT) );

  memcpy(outputClusterProminence, clusterProminenceImage->GetBufferPointer(),
         clusterProminenceImage->GetBufferedRegion().GetNumberOfPixels()*sizeof(OT) );
 
}

//-----------------------------------------------------------------------------
vtkITKHaralickTextureFeaturesImageFilter::vtkITKHaralickTextureFeaturesImageFilter()
{
   this->SetRadius(3);
   this->SetNumberOfBins(128);
   this->SetTextureType(0); // 0 - fine, 1 - medium, 2 - coarse
   this->SetNormalizeTextures(1);
   this->SetApplyMedianFilter(1);
   this->SetNormMaximum(256);
   this->SetNormMinimum(0);
 
   // this->SetNumberOfInputPorts(1);
   this->SetNumberOfOutputPorts(7); // as we will have 7 different outputs for texture
}

//-----------------------------------------------------------------------------
template< class IT1>
void ExecuteHaralickTextureFeatures( vtkITKHaralickTextureFeaturesImageFilter *self,
          vtkImageData *input1,
	  vtkImageData *input2,
	  vtkImageData *outData1,
	  vtkImageData *outData2,
	  vtkImageData *outData3,
	  vtkImageData *outData4,
	  vtkImageData *outData5,
	  vtkImageData *outData6,
	  vtkImageData *outData7,
          IT1 *)
{
  int outExt[6];
  int dims[3];
  double spacing[3], origin[3];

  input1->GetDimensions(dims);
  input1->GetOrigin(origin);
  input1->GetSpacing(spacing);
  input1->GetExtent(outExt);

  void *inPtr1 = input1->GetScalarPointerForExtent(outExt);
  void *inPtr2 =  input2->GetScalarPointerForExtent(outExt);

  input1->GetWholeExtent(outExt);
  outData1->SetExtent(outExt);
  outData1->SetOrigin(origin);
  outData1->SetSpacing(spacing);
  outData1->SetDimensions(dims);
  outData1->AllocateScalars();

  outData2->SetExtent(outExt);
  outData2->SetOrigin(origin);
  outData2->SetSpacing(spacing);
  outData2->SetDimensions(dims);
  outData2->AllocateScalars();

  outData3->SetExtent(outExt);
  outData3->SetOrigin(origin);
  outData3->SetSpacing(spacing);
  outData3->SetDimensions(dims);
  outData3->AllocateScalars();

  outData4->SetExtent(outExt);
  outData4->SetOrigin(origin);
  outData4->SetSpacing(spacing);
  outData4->SetDimensions(dims);
  outData4->AllocateScalars();

  outData5->SetExtent(outExt);
  outData5->SetOrigin(origin);
  outData5->SetSpacing(spacing);
  outData5->SetDimensions(dims);
  outData5->AllocateScalars();

  outData6->SetExtent(outExt);
  outData6->SetOrigin(origin);
  outData6->SetSpacing(spacing);
  outData6->SetDimensions(dims);
  outData6->AllocateScalars();
  
  outData7->SetExtent(outExt);
  outData7->SetOrigin(origin);
  outData7->SetSpacing(spacing);
  outData7->SetDimensions(dims);
  outData7->AllocateScalars();

  input1->GetExtent(outExt);

  void *outPtr1 = outData1->GetScalarPointerForExtent(outExt);
  void *outPtr2 = outData2->GetScalarPointerForExtent(outExt);
  void *outPtr3 = outData3->GetScalarPointerForExtent(outExt);
  void *outPtr4 = outData4->GetScalarPointerForExtent(outExt);
  void *outPtr5 = outData5->GetScalarPointerForExtent(outExt);
  void *outPtr6 = outData6->GetScalarPointerForExtent(outExt);
  void *outPtr7 = outData7->GetScalarPointerForExtent(outExt);
  
  // set up progress callback
  //itk::CStyleCommand::Pointer progressCommand = itk::CStyleCommand::New();
  //progressCommand->SetClientData(static_cast<void *>(self));
  //progressCommand->SetCallback(vtkITKImageGrowCutHandleProgressEvent );


  if(self->TextureType == 0)
  {
     self->SetRadius(3);
     self->SetNumberOfBins(128);
  }
  else if(self->TextureType == 1)
  {
    self->SetRadius(5);
    self->SetNumberOfBins(128);
  }
  else
  {
     self->SetRadius(9);
     self->SetNumberOfBins(128);
  }

  std::cout<<" In Execute Main "<<std::endl;

  if((input2->GetScalarType() != VTK_UNSIGNED_SHORT) ||
     (input2->GetScalarType() != VTK_UNSIGNED_CHAR) ||
     (input2->GetScalarType() != VTK_UNSIGNED_LONG) ||
     (input2->GetScalarType() != VTK_SHORT) ||
     (input2->GetScalarType() != VTK_CHAR) ||
     (input2->GetScalarType() != VTK_LONG) )
    {

        vtkImageCast *imageCaster1 = vtkImageCast::New();
        imageCaster1->SetInput(input2);
        imageCaster1->SetOutputScalarTypeToShort();

	vtkImageCast *imageCaster = vtkImageCast::New();
	vtkImageCast *imageCasterEnergy = vtkImageCast::New();
	vtkImageCast *imageCasterEntropy = vtkImageCast::New();
	vtkImageCast *imageCasterCorr = vtkImageCast::New();
	vtkImageCast *imageCasterDiffMoment = vtkImageCast::New();
	vtkImageCast *imageCasterInertia = vtkImageCast::New();
	vtkImageCast *imageCasterCShade = vtkImageCast::New();
	vtkImageCast *imageCasterCProminence = vtkImageCast::New();
	
	if(input1->GetScalarType() != VTK_SHORT)
	{
	  imageCaster->SetInput(input1);
	  imageCaster->SetOutputScalarTypeToShort();
	  
	  imageCasterEnergy->SetInput(outData1);
	  imageCasterEnergy->SetOutputScalarTypeToShort();

	  imageCasterEntropy->SetInput(outData2);
	  imageCasterEntropy->SetOutputScalarTypeToShort();

	  imageCasterCorr->SetInput(outData3);
	  imageCasterCorr->SetOutputScalarTypeToShort();

	  imageCasterDiffMoment->SetInput(outData4);
	  imageCasterDiffMoment->SetOutputScalarTypeToShort();

	  imageCasterInertia->SetInput(outData5);
	  imageCasterInertia->SetOutputScalarTypeToShort();

	  imageCasterCShade->SetInput(outData6);
	  imageCasterCShade->SetOutputScalarTypeToShort();

	  imageCasterCProminence->SetInput(outData7);
	  imageCasterCProminence->SetOutputScalarTypeToShort();
        }

        vtkITKImageHaralickTextureFeaturesExecute3D(input1, (short*)(inPtr1), 
		   (short*)(inPtr2), (short*)(outPtr1), 
		   (short*)(outPtr2), (short*)(outPtr3),
		   (short*)(outPtr4), (short*)(outPtr5),
		   (short*)(outPtr6), (short*)(outPtr7),
		   self->Radius, self->NumberOfBins,
		   (self->NormalizeTextures == 1 ? true : false),
		   (self->ApplyMedianFilter == 1 ? true : false),
		   self->NormMaximum, 
		   self->NormMinimum, self->m_IndexLabelMap, self->m_MeanEnergy, 
		   self->m_MeanEntropy, self->m_MeanCorrelation, self->m_MeanInertia,
                   self->m_MeanDifferenceMoment, self->m_MeanClusterShade,
		   self->m_MeanClusterProminence, self->m_StdEnergy,
                   self->m_StdEntropy, self->m_StdCorrelation,
                   self->m_StdDifferenceMoment, self->m_StdInertia,
                   self->m_StdClusterShade, self->m_StdClusterProminence,
                   self->m_KurtosisEnergy, self->m_KurtosisEntropy,
		   self->m_KurtosisCorrelation, self->m_KurtosisDifferenceMoment,
                   self->m_KurtosisInertia, self->m_KurtosisClusterShade,
                   self->m_KurtosisClusterProminence, self->m_SkewEnergy, 
		   self->m_SkewEntropy, self->m_SkewCorrelation, self->m_SkewDifferenceMoment, 
		   self->m_SkewInertia, self->m_SkewClusterShade, self->m_SkewClusterProminence, 
		   self->m_MedianEnergy, self->m_MedianEntropy, self->m_MedianCorrelation, 
		   self->m_MedianDifferenceMoment, self->m_MedianInertia, self->m_MedianClusterShade, 
		   self->m_MedianClusterProminence);

	self->NumberLabels = self->m_IndexLabelMap.size();
	std::cout<<"number of labels "<<self->NumberLabels<<std::endl;
        imageCaster1->Delete();

	imageCaster->Delete();
	imageCasterEnergy->Delete();
	imageCasterEntropy->Delete();
	imageCasterCorr->Delete();
	imageCasterDiffMoment->Delete();
	imageCasterInertia->Delete();
	imageCasterCShade->Delete();
	imageCasterCProminence->Delete();
    }
}

//-----------------------------------------------------------------------------
void vtkITKHaralickTextureFeaturesImageFilter::ExecuteData(
 vtkDataObject *outData)
{
  vtkImageData *input1 = GetInput(0);
  vtkImageData *input2 = GetInput(1);

  vtkImageData * out1 = vtkImageData::SafeDownCast(outData);

  // get the remaining outputs
  vtkImageData *out2 = vtkImageData::SafeDownCast(this->GetOutputDataObject(1)); //this->m_EntropyImage; //this->GetOutput(2);
  vtkImageData *out3 = vtkImageData::SafeDownCast(this->GetOutputDataObject(2)); //this->m_CorrelationImage; //this->GetOutput(3);
  vtkImageData *out4 = vtkImageData::SafeDownCast(this->GetOutputDataObject(3)); //this->m_DifferenceMomentImage; //this->GetOutput(4);
  vtkImageData *out5 = vtkImageData::SafeDownCast(this->GetOutputDataObject(4)); //this->m_InertiaImage; //this->GetOutput(5);
  vtkImageData *out6 = vtkImageData::SafeDownCast(this->GetOutputDataObject(5)); //this->m_ClusterShadeImage; //this->GetOutput(6);
  vtkImageData *out7 = vtkImageData::SafeDownCast(this->GetOutputDataObject(6)); //this->m_ClusterProminenceImage; //this->GetOutput(7);

  switch(input1->GetScalarType() ) {
    vtkTemplateMacro( ExecuteHaralickTextureFeatures(this, input1, input2,					     out1, out2, out3, out4, out5, out6, out7,
             static_cast< VTK_TT*>(0)));
    break;
  }
}

//-----------------------------------------------------------------------------
void vtkITKHaralickTextureFeaturesImageFilter::ExecuteInformation()
{
  this->Superclass::ExecuteInformation();
}

//-----------------------------------------------------------------------------
void vtkITKHaralickTextureFeaturesImageFilter::ExecuteInformation(
  vtkImageData **inputs, vtkImageData **outputs)
{
  std::cout<<" In Haralick Texture Filter Execute Information..."<<std::endl;
  vtkImageData* output;
  if (inputs[1] == NULL)
  {
    std::cout<<" second input is not set "<<std::endl;
    return;
  }

  for (int i = 0; i < this->NumberOfOutputs; ++i)
  {
    output = vtkImageData::SafeDownCast(this->GetOutputDataObject(i));
    if(output)
    {
      output->CopyTypeSpecificInformation(inputs[0]);
    }
  }
}

//-----------------------------------------------------------------------------
void vtkITKHaralickTextureFeaturesImageFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "Radius : " << this->Radius << std::endl;
}
