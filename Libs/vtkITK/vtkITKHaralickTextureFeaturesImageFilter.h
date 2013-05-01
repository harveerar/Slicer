/*=================================================================================================================
 Copyright Memorial Sloan-Kettering Cancer Center (MSKCC) All Rights Reserved.
 
 Module: HaralickTextureFeaturesImageFilterWrapper
 Author: Harini Veeraraghavan
 Email: veerarah at mskcc dot org
 Version: $Revision: 1.2$
 Date: March 22nd 2013  

  ==================================================================================================================*/




#ifndef __vtkITKHaralickTextureFeaturesImageFilter_h
#define __vtkITKHaralickTextureFeaturesImageFilter_h

#include "vtkITK.h"
#include "vtkPolyData.h"
#include <vector>

// VTK includes
#include <vtkImageMultipleInputOutputFilter.h>

class vtkImageData;

typedef enum{
  Energy,
  Entropy,
  Correlation,
  DifferenceMoment,
  Interia,
  ClusterShade,
  ClusterProminence
}HaralickFeatures;

typedef enum{
  fine,
  medium,
  coarse
}ParameterSetting;

/// \brief- Wrapper class around itk::HaralickTextureFeaturesImageFilter
///
/// Computes Morphological textures in a given ROI specified by a label map
/// automatically computes the region of interest from the label map and computes the Haralick 
/// textures per voxel inside the label map
///
/// Usage: SetInput1 is the input feature/intensity image (required)
/// SetInput2 takes the label image (the label mask specifying ROI) image (required)
///
/// Produces multiple outputs 1 for each texture feature
/// GetOutput1 - energy
/// GetOutput2 - entropy
/// GetOutput3 - Correlation
/// GetOutput4 - DifferenceMoment
/// GetOutput5 - Intertia
/// GetOutput6 - ClusterShade
/// GetOutput7 - ClusterProminence
///
/// This filter is implemented only for scalar images gray scale images. 
/// The current implementation supports n-class segmentation.
class VTK_ITK_EXPORT vtkITKHaralickTextureFeaturesImageFilter : public vtkImageMultipleInputOutputFilter 
{
public:

  static vtkITKHaralickTextureFeaturesImageFilter *New();
  vtkTypeRevisionMacro(vtkITKHaralickTextureFeaturesImageFilter,vtkImageMultipleInputOutputFilter );
  void PrintSelf(ostream& os, vtkIndent indent);

  /// Methods to set/get texture radius
  vtkSetMacro( Radius, double);
  vtkGetMacro( Radius, double);

  vtkSetMacro( NumberOfBins, double);
  vtkGetMacro( NumberOfBins, double);

  vtkSetMacro( TextureType, int);
  vtkGetMacro( TextureType, int);

  vtkSetMacro(NormalizeTextures, int);
  vtkGetMacro(NormalizeTextures, int);
  
  vtkSetMacro(ApplyMedianFilter, int);
  vtkGetMacro(ApplyMedianFilter, int);

  vtkSetMacro(NormMaximum, float);
  vtkGetMacro(NormMaximum, float);

  vtkSetMacro(NormMinimum, float);
  vtkGetMacro(NormMinimum, float);

  vtkSetMacro(NumberLabels, int);
  vtkGetMacro(NumberLabels, int);
  
  // outputs for various texture statistics
  float GetMeanEnergyForLabelIndex(int index)
  {
    if (index < m_MeanEnergy.size())
    { 
      return m_MeanEnergy[index];
    }
    return 0.0;
  }
  float GetMeanEntropyForLabelIndex(int index)
  {
    if(index < m_MeanEntropy.size())
    { 
      return m_MeanEntropy[index];
    }
    return 0.0;
  }

  float GetMeanCorrelationForLabelIndex(int index)
  {
    if (index < m_MeanCorrelation.size())
    {
      return m_MeanCorrelation[index];
    }
    return 0.0;
  }
  float GetMeanDifferenceMomentForLabelIndex(int index)
  {
    if(index < m_MeanDifferenceMoment.size())
    { 
      return m_MeanDifferenceMoment[index];
    }
    return 0.0;
  }
  float GetMeanInertiaForLabelIndex(int index)
  {
    if(index < m_MeanInertia.size())
    { 
      return m_MeanInertia[index];
    }
    return 0.0;
  }
  float GetMeanClusterShadeForLabelIndex(int index)
  {
    if(index < m_MeanClusterShade.size())
    { 
      return m_MeanClusterShade[index];
    }
    return 0.0;
  }
  float GetMeanClusterProminenceForLabelIndex(int index)
  {
    if(index < m_MeanClusterProminence.size())
    { 
      return m_MeanClusterProminence[index];
    }
    return 0.0;
  }
  
  float GetKurtosisEnergyForLabelIndex(int index)
  {
    if(index < m_KurtosisEnergy.size())
    { 
      return m_KurtosisEnergy[index];
    }
    return 0.0;
  }
  float GetKurtosisEntropyForLabelIndex(int index)
  {
    if(index < m_KurtosisEntropy.size())
    { 
      return m_KurtosisEntropy[index];
    }
    return 0.0;
  }
  float GetKurtosisCorrelationForLabelIndex(int index)
  {
    if(index < m_KurtosisCorrelation.size())
    {
      return m_KurtosisCorrelation[index];
    }
    return 0.0;
  }
  float GetKurtosisInertiaForLabelIndex(int index)
  {
    if(index < m_KurtosisInertia.size())
    { 
      return m_KurtosisInertia[index];
    }
    return 0.0;
  }
  float GetKurtosisDifferenceMomentForLabelIndex(int index)
  {
    if(index < m_KurtosisDifferenceMoment.size())
    { 
      return m_KurtosisDifferenceMoment[index];
    }
    return 0.0;
  }
  float GetKurtosisClusterShadeForLabelIndex(int index)
  {
    if(index < m_KurtosisClusterShade.size())
    { 
      return m_KurtosisClusterShade[index];
    }
    return 0.0;
  }
  float GetKurtosisClusterProminenceForLabelIndex(int index)
  {
    if(index < m_KurtosisClusterProminence.size())
    { 
      return m_KurtosisClusterProminence[index];
    }
    return 0.0;
  }
  
  float GetStdEnergyForLabelIndex(int index)
  {
    if(index < m_StdEnergy.size())
    { 
      return m_StdEnergy[index];
    }
    return 0.0;
  }

  float GetSkewEnergyForLabelIndex(int index)
  {
    if(index < m_SkewEnergy.size())
    { 
      return m_SkewEnergy[index];
    }
    return 0.0;
  }

  float GetMedianEnergyForLabelIndex(int index)
  {
    if(index < m_MedianEnergy.size())
    { 
      return m_MedianEnergy[index];
    }
    return 0.0;
  }

  float GetStdEntropyForLabelIndex(int index)
  {
    if(index < m_StdEntropy.size())
    { 
      return m_StdEntropy[index];
    }
    return 0.0;
  }

  float GetSkewEntropyForLabelIndex(int index)
  {
    if(index < m_SkewEntropy.size())
    { 
      return m_SkewEntropy[index];
    }
    return 0.0;
  }

  float GetMedianEntropyForLabelIndex(int index)
  {
    if(index < m_MedianEntropy.size())
    { 
      return m_MedianEntropy[index];
    }
    return 0.0;
  }

  float GetStdCorrelationForLabelIndex(int index)
  {
    if(index < m_StdCorrelation.size())
    { 
      return m_StdCorrelation[index];
    }
    return 0.0;
  }

  float GetSkewCorrelationForLabelIndex(int index)
  {
    if(index < m_SkewCorrelation.size())
    { 
      return m_SkewCorrelation[index];
    }
    return 0.0;
  }

  float GetMedianCorrelationForLabelIndex(int index)
  {
    if(index < m_MedianCorrelation.size())
    { 
      return m_MedianCorrelation[index];
    }
    return 0.0;
  }

  float GetStdInertiaForLabelIndex(int index)
  {
    if(index < m_StdInertia.size())
    { 
      return m_StdInertia[index];
    }
    return 0.0;
  }

  float GetSkewInertiaForLabelIndex(int index)
  {
    if(index < m_SkewInertia.size())
    { 
      return m_SkewInertia[index];
    }
    return 0.0;
  }

  float GetMedianInertiaForLabelIndex(int index)
  {
    if(index < m_MedianInertia.size())
    { 
      return m_MedianInertia[index];
    }
    return 0.0;
  }

  float GetStdDifferenceMomentForLabelIndex(int index)
  {
    if(index < m_StdDifferenceMoment.size())
    { 
      return m_StdDifferenceMoment[index];
    }
    return 0.0;
  }

  float GetSkewDifferenceMomentForLabelIndex(int index)
  {
    if(index < m_SkewDifferenceMoment.size())
    { 
      return m_SkewDifferenceMoment[index];
    }
    return 0.0;
  }

  float GetMedianDifferenceMomentForLabelIndex(int index)
  {
    if(index < m_MedianDifferenceMoment.size())
    { 
      return m_MedianDifferenceMoment[index];
    }
    return 0.0;
  }

  float GetStdClusterShadeForLabelIndex(int index)
  {
    if(index < m_StdClusterShade.size())
    { 
      return m_StdClusterShade[index];
    }
    return 0.0;
  }

  float GetSkewClusterShadeForLabelIndex(int index)
  {
    if(index < m_SkewClusterShade.size())
    { 
      return m_SkewClusterShade[index];
    }
    return 0.0;
  }

  float GetMedianClusterShadeForLabelIndex(int index)
  {
    if(index < m_MedianClusterShade.size())
    { 
      return m_MedianClusterShade[index];
    }
    return 0.0;
  }

  float GetStdClusterProminenceForLabelIndex(int index)
  {
    if(index < m_StdClusterProminence.size())
    { 
      return m_StdClusterProminence[index];
    }
      return 0.0;
  }

  float GetSkewClusterProminenceForLabelIndex(int index)
  {
    if(index < m_SkewClusterProminence.size())
    { 
      return m_SkewClusterProminence[index];
    }
      return 0.0;
  }

  float GetMedianClusterProminenceForLabelIndex(int index)
  {
    if(index < m_MedianClusterProminence.size())
    { 
      return m_MedianClusterProminence[index];
    }
      return 0.0;
  }

  int GetLabelForIndex(int index)
  {
    if(index < m_IndexLabelMap.size())
    { 
      return m_IndexLabelMap[index];
    }
    return 0.0;
  }

  // return histogram of each feature
  vtkPolyData* GetEnergyHistogramForIndex( int index)
  {
    if (index < m_EnergyHistogram.size())
      return m_EnergyHistogram[index];
  }
  
  vtkPolyData* GetEntropyHistogramForIndex( int index)
  {
    if (index < m_EntropyHistogram.size())
      return m_EntropyHistogram[index];
  }

  vtkPolyData* GetCorrelationHistogramForIndex(int index)
  {
    if(index < m_CorrelationHistogram.size())
      return m_CorrelationHistogram[index];
  }

  vtkPolyData* GetDifferenceMomentHistogramForIndex(int index)
  {
    if(index < m_DifferenceMomentHistogram.size())
      return m_DifferenceMomentHistogram[index];
  }

  vtkPolyData* GetInertiaHistogramForIndex(int index)
  {
    if(index < m_InertiaHistogram.size())
      return m_InertiaHistogram[index];
  }

  /*  vtkPolyData* GetClusterShadeHistogram()
  {
    if(index < m_ClusterShadeHistogram.size())
      return m_ClusterShadeHistogram[index];
  }

  vtkPolyData* GetClusterProminenceHistogram()
  {
    if(index < m_ClusterProminenceHistogram.size())
      return m_ClusterProminenceHistogram[index];
      }*/

  
  
  // For future -- return a 5-bin histogram for each feature

public:
  double     Radius;
  double     NumberOfBins;
  int        TextureType;   // 0 -fine 1 - medium, 2 - coarse
  int        NormalizeTextures;  // 1 for yes and 0 for no

  int        ApplyMedianFilter;

  float      NormMaximum;

  float      NormMinimum;

  int        NumberLabels;

  std::vector< int >    m_IndexLabelMap;
  std::vector< float >  m_MeanEnergy;
  std::vector< float >  m_MeanEntropy;
  std::vector< float >  m_MeanCorrelation;
  std::vector< float >  m_MeanDifferenceMoment;
  std::vector< float >  m_MeanInertia;
  std::vector< float >  m_MeanClusterShade;
  std::vector< float >  m_MeanClusterProminence;

  std::vector< float >  m_MedianEnergy;
  std::vector< float >  m_MedianEntropy;
  std::vector< float >  m_MedianCorrelation;
  std::vector< float >  m_MedianDifferenceMoment;
  std::vector< float >  m_MedianInertia;
  std::vector< float >  m_MedianClusterShade;
  std::vector< float >  m_MedianClusterProminence;

  std::vector< float >  m_StdEnergy;
  std::vector< float >  m_StdEntropy;
  std::vector< float >  m_StdCorrelation;
  std::vector< float >  m_StdDifferenceMoment;
  std::vector< float >  m_StdInertia;
  std::vector< float >  m_StdClusterShade;
  std::vector< float >  m_StdClusterProminence;

  std::vector< float >  m_KurtosisEnergy;
  std::vector< float >  m_KurtosisEntropy;
  std::vector< float >  m_KurtosisCorrelation;
  std::vector< float >  m_KurtosisDifferenceMoment;
  std::vector< float >  m_KurtosisInertia;
  std::vector< float >  m_KurtosisClusterShade;
  std::vector< float >  m_KurtosisClusterProminence;

  std::vector< float >  m_SkewEnergy;
  std::vector< float >  m_SkewEntropy;
  std::vector< float >  m_SkewCorrelation;
  std::vector< float >  m_SkewDifferenceMoment;
  std::vector< float >  m_SkewInertia;
  std::vector< float >  m_SkewClusterShade;
  std::vector< float >  m_SkewClusterProminence;


  std::vector< vtkPolyData* > m_EnergyHistogram;
  std::vector< vtkPolyData* > m_EntropyHistogram;
  std::vector< vtkPolyData* > m_CorrelationHistogram;
  std::vector< vtkPolyData* > m_DifferenceMomentHistogram;
  std::vector< vtkPolyData* > m_InertiaHistogram;
  std::vector< vtkPolyData* > m_ClusterShadeHistogram;
  std::vector< vtkPolyData* > m_ClusterProminenceHistogram;

  std::vector< vtkPolyData* > m_EnergyFeaturePoints;
  std::vector< vtkPolyData* > m_EntropyFeaturePoints;
  std::vector< vtkPolyData* > m_CorrelationFeaturePoints;
  std::vector< vtkPolyData* > m_DifferenceMomentFeaturePoints;
  std::vector< vtkPolyData* > m_InertiaFeaturePoints;
  std::vector< vtkPolyData* >  m_ClusterShadeFeaturePoints;
  std::vector< vtkPolyData* >  m_ClusterProminenceFeaturePoints;

 protected:
  vtkITKHaralickTextureFeaturesImageFilter();
  ~vtkITKHaralickTextureFeaturesImageFilter(){}

  virtual void ExecuteData(vtkDataObject *outData);

  /// Override ExecuteInformation so that the second input is used to
  /// define the output information (input gestures and output
  /// segmentation images should be same image type)
  virtual void ExecuteInformation(vtkImageData **, vtkImageData **);

  /// Need to provide ExecuteInformation() or it will be hidden by the
  /// override to ExecuteInformation(vtkImageData**, vtkImageData**)
  virtual void ExecuteInformation();

private:
  vtkITKHaralickTextureFeaturesImageFilter(const vtkITKHaralickTextureFeaturesImageFilter&);  // Not implemented.
  void operator=(const vtkITKHaralickTextureFeaturesImageFilter&);  // Not implemented.

};

#endif
