/*
 * stkITKTools.hxx
 *
 *  Created on: 6 Aug 2015
 *      Author: lewish
 */

#ifndef __stkITKTools_hxx
#define __stkITKTools_hxx


#include  "stkITKTools.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkDataObject.h"
#include "itkImageFileWriter.h"
#include "itkRegularStepGradientDescentOptimizer.h"
namespace stk{

ITKTools::ITKTools(){}
ITKTools::~ITKTools(){}

typedef itk::Image<float, 2> ImageType;
ImageType::Pointer ITKTools::GaussianBlur(ImageType::Pointer itkInput){

	typedef itk::Image<float, 2> ImageType;
	typedef itk::DiscreteGaussianImageFilter< ImageType, ImageType> FilterType;

	FilterType::Pointer filter = FilterType::New();
	filter->SetInput(itkInput);


	filter->SetVariance( 15 );
	//filter->SetMaximumKernelWidth( 20 );
	filter->Update();
	return filter->GetOutput();




}



ImageType::Pointer ITKTools::EdgeDetection(ImageType::Pointer itkInput){


	typedef itk::SobelEdgeDetectionImageFilter <ImageType, ImageType>  SobelFilterType;
	SobelFilterType::Pointer sobelFilter = SobelFilterType::New();
	sobelFilter->SetInput(itkInput);
	sobelFilter->Update();

	return sobelFilter->GetOutput();




}

typedef itk::TranslationTransform< double, 2 > TransformType;
ImageType::Pointer ITKTools::ImageRegistration(ImageType::Pointer fixedImage, ImageType::Pointer movingImage){


	typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
	typedef itk::MeanSquaresImageToImageMetric< ImageType, ImageType >    MetricType;
	typedef itk:: LinearInterpolateImageFunction< ImageType, double >    InterpolatorType;
	typedef itk::ImageRegistrationMethod< ImageType, ImageType >    RegistrationType;


	MetricType::Pointer         metric        = MetricType::New();
	TransformType::Pointer      transform     = TransformType::New();
	OptimizerType::Pointer      optimizer     = OptimizerType::New();
	InterpolatorType::Pointer   interpolator  = InterpolatorType::New();
	RegistrationType::Pointer   registration  = RegistrationType::New();


	registration->SetMetric(        metric        );
	registration->SetOptimizer(     optimizer     );
	registration->SetTransform(     transform     );
	registration->SetInterpolator(  interpolator  );


	registration->SetFixedImage(fixedImage);
	registration->SetMovingImage(movingImage);

	registration->SetFixedImageRegion(fixedImage->GetLargestPossibleRegion() );
	typedef RegistrationType::ParametersType ParametersType;
	ParametersType initialParameters( transform->GetNumberOfParameters() );

	initialParameters[0] = 0.0;  // Initial offset along X
	initialParameters[1] = 0.0;  // Initial offset along Y

	registration->SetInitialTransformParameters( initialParameters );

	optimizer->SetMaximumStepLength( 7.5);
	optimizer->SetMinimumStepLength( 0.01 );

	// Set a stopping criterion
	optimizer->SetNumberOfIterations( 200 );

	registration->Update();

	ParametersType finalParameters = registration->GetLastTransformParameters();



	const double TranslationAlongX = finalParameters[0];
	const double TranslationAlongY = finalParameters[1];



	const unsigned int numberOfIterations = optimizer->GetCurrentIteration();



	const double bestValue = optimizer->GetValue();

	// Print out results
	//
	std::cout << "Result = " << std::endl;
	std::cout << " Translation X = " << TranslationAlongX  << std::endl;
	std::cout << " Translation Y = " << TranslationAlongY  << std::endl;
	std::cout << " Iterations    = " << numberOfIterations << std::endl;
	std::cout << " Metric value  = " << bestValue          << std::endl;


	typedef itk::ResampleImageFilter< ImageType, ImageType >    ResampleFilterType;

	//  A resampling filter is created and the moving image is connected as  its input.

	ResampleFilterType::Pointer resampler = ResampleFilterType::New();
	resampler->SetInput( movingImage);
	resampler->SetTransform( registration->GetOutput()->Get() );
	resampler->SetSize( fixedImage->GetLargestPossibleRegion().GetSize() );
	resampler->SetOutputOrigin(  fixedImage->GetOrigin() );
	resampler->SetOutputSpacing( fixedImage->GetSpacing() );
	resampler->SetOutputDirection( fixedImage->GetDirection() );
	resampler->SetDefaultPixelValue( 0.2 );

	resampler->Update();

	return resampler->GetOutput();
}


ImageType::Pointer ITKTools::ImageRotation(ImageType::Pointer fixedImage, ImageType::Pointer movingImage){

	typedef itk::AffineTransform< double, 2 > TransformType;
	typedef itk::RegularStepGradientDescentOptimizer       OptimizerType;
	typedef itk::MeanSquaresImageToImageMetric<ImageType, ImageType >    MetricType;
	typedef itk:: LinearInterpolateImageFunction<ImageType,	double >    InterpolatorType;
	typedef itk::ImageRegistrationMethod< ImageType, ImageType >    RegistrationType;

	// Create components
	MetricType::Pointer         metric        = MetricType::New();
	TransformType::Pointer      transform     = TransformType::New();
	OptimizerType::Pointer      optimizer     = OptimizerType::New();
	InterpolatorType::Pointer   interpolator  = InterpolatorType::New();
	RegistrationType::Pointer   registration  = RegistrationType::New();

	// Each component is now connected to the instance of the registration method.
	registration->SetMetric(        metric        );
	registration->SetOptimizer(     optimizer     );
	registration->SetTransform(     transform     );
	registration->SetInterpolator(  interpolator  );
	registration->SetFixedImage(fixedImage);
	registration->SetMovingImage(movingImage);

	registration->SetFixedImageRegion(
			fixedImage->GetLargestPossibleRegion() );

	//  Initialize the transform
	typedef RegistrationType::ParametersType ParametersType;
	ParametersType initialParameters( transform->GetNumberOfParameters() );

	// rotation matrix
	initialParameters[0] = 1;  // R(0,0)
	initialParameters[1] = 0;  // R(0,1)
	initialParameters[2] = 0;  // R(1,0)
	initialParameters[3] = 1;  // R(1,1)

	// translation vector
	initialParameters[4] = 0.0;
	initialParameters[5] = 0.0;

	registration->SetInitialTransformParameters( initialParameters );

	optimizer->SetMaximumStepLength( 0.1 ); // If this is set too high, you will get a
	//"itk::ERROR: MeanSquaresImageToImageMetric(0xa27ce70): Too many samples map outside moving image buffer: 1818 / 10000" error

	optimizer->SetMinimumStepLength( 0.001 );

	// Set a stopping criterion
	optimizer->SetNumberOfIterations( 200 );
	registration->Update();

	ParametersType finalParameters = registration->GetLastTransformParameters();
	std::cout << "Final parameters: " << finalParameters << std::endl;

	//  The value of the image metric corresponding to the last set of parameters
	//  can be obtained with the \code{GetValue()} method of the optimizer.

	const double bestValue = optimizer->GetValue();

	// Print out results
	//
	std::cout << "Result = " << std::endl;
	std::cout << " Metric value  = " << bestValue          << std::endl;

	//  It is common, as the last step of a registration task, to use the
	//  resulting transform to map the moving image into the fixed image space.
	//  This is easily done with the \doxygen{ResampleImageFilter}.

	typedef itk::ResampleImageFilter<ImageType, ImageType >    ResampleFilterType;

	ResampleFilterType::Pointer resampler = ResampleFilterType::New();
	resampler->SetInput( movingImage);



	resampler->SetTransform( registration->GetOutput()->Get() );


	resampler->SetSize( fixedImage->GetLargestPossibleRegion().GetSize() );
	resampler->SetOutputOrigin(  fixedImage->GetOrigin() );
	resampler->SetOutputSpacing( fixedImage->GetSpacing() );
	resampler->SetOutputDirection( fixedImage->GetDirection() );
	resampler->SetDefaultPixelValue( 100 );

	resampler->Update();

	return resampler->GetOutput();

}
}





#endif /* MODULES_IMAGECALCULATIONS_ITKTOOLS_INCLUDE_STKITKTOOLS_HXX_ */
