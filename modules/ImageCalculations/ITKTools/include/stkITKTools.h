/*
 * stkITKTools.h
 *
 *  Created on: 6 Aug 2015
 *      Author: lewish
 */

#ifndef __stkITKTools_h
#define __stkITKTools_h

#include <memory>
#include <iterator>
#include <vector>

#include "stkImage.h"
#include "stkImageStack.h"
#include "itkImage.h"
#include "itkImageRegionIterator.h"
#include "itkImageFileWriter.h"
#include "itkImageBase.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkDataObject.h"

#include "itkImageFileWriter.h"
#include "itkImageToImageFilter.h"
#include "itkDataObject.h"

#include "itkImageSource.h"
#include "itkSobelEdgeDetectionImageFilter.h"
#include "itkImageRegistrationMethod.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkInterpolateImageFunction.h"
#include "itkTranslationTransform.h"


#include "itkCastImageFilter.h"
#include "itkEllipseSpatialObject.h"
#include "itkImage.h"
#include "itkImageRegistrationMethod.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkResampleImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkSpatialObjectToImageFilter.h"
#include "itkTranslationTransform.h"
#include "itkAffineTransform.h"

namespace stk{


class ITKTools{

public:

	ITKTools();
	virtual ~ITKTools();


typedef itk::Image<float, 2> ImageType;
ImageType::Pointer GaussianBlur(ImageType::Pointer itkInput);



ImageType::Pointer EdgeDetection(ImageType::Pointer itkInput);

typedef itk::TranslationTransform< double, 2 > TransformType;
ImageType::Pointer ImageRegistration(ImageType::Pointer fixedImage, ImageType::Pointer movingImage);


ImageType::Pointer ImageRotation(ImageType::Pointer fixedImage, ImageType::Pointer movingImage);


};




}
#include "stkITKTools.hxx"

#endif /* MODULES_IMAGECALCULATIONS_ITKTOOLS_INCLUDE_STKITKTOOLS_H_ */
