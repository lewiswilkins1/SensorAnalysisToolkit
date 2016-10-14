/*
 * EdgeDetector.cxx
 *
 *  Created on: 3 Sep 2015
 *      Author: lewish
 */




#include <iostream>
#include <sstream>


#include <chrono>
#include <ctime>
#include <ratio>
#include <cmath>


#include "stkImageHistogram.h"//write to histogram
#include "stkRawImageIO.h"//write pedestal in the frames
#include "stkRawImageStackIO.h"//load in the frames
#include "stkImageStack.h"//hold the files to be loaded
#include "stkImage.h"//will hold the result
#include "stkImageSum.h"//used to sum the stack
#include "stkImageDivision.h"//used to divide through to get final result
#include "stkImageVariance.h"//used for variance
#include "stkImageMask.h"
#include "stkImageBadPixelAverage.h"
#include "stkImageMinus.h"
#include "stkImageResize.h"
#include "stkITKAdapter.h"
#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "stkITKTools.h"
#include "itkImageFileWriter.h"

int main(int argc, const char** argv){

	std::string pedfilePath,fileNameAndFormat, lightfilePath, lightfileNameAndFormat, filePathandfilename;
	int startingframe, numOfFrames, rows, cols, framesize;
	std::stringstream inputs;
	for(int iArg=1; iArg < argc; iArg++)
	{
		inputs << argv[iArg] << ' ';//get the arguments
	}

	inputs >>  filePathandfilename >>   rows >> cols;//write the inputs
	framesize = cols*rows;




	stk::IO<float> myIO;
	std::shared_ptr< stk::Image<float> > myImage( new stk::Image<float> );
	myImage->Initialise( myIO.ReadImage( filePathandfilename ), rows, cols );


	std::shared_ptr<stk::Image<float> > myResultResize (new stk::Image<float>(1024,1024));

	stk::ImageResize mySize;
	mySize.ResizeImage(myImage , myResultResize);

	std::cout<<myResultResize->GetPixelAt(3874);

	std::shared_ptr<stk::Image<float> > myResultResize2 (new stk::Image<float>(1024,1024));



		typedef itk::Image< float, 2 >  ImageType;
		ImageType::Pointer itkImage = ImageType::New();
		ImageType::Pointer itkBlurredImage = ImageType::New();
		ImageType::Pointer itkSobelImage = ImageType::New();



		ImageType::RegionType region;
		ImageType::IndexType start;


		ImageType::SizeType size;
		start[0] = 0;
		start[1] = 0;
		size[0] = rows;
		size[1] = cols;
		region.SetSize(size);
		region.SetIndex(start);
		itkBlurredImage->SetRegions(region);
		itkBlurredImage->Allocate();

		itkSobelImage->SetRegions(region);
		itkSobelImage->Allocate();






		stk::ITKAdapter adapter;
		adapter.STKtoITK(itkImage, myResultResize, 1024, 1024);


		stk::ITKTools gaussian;

		itkBlurredImage = gaussian.GaussianBlur(itkImage);
		itkBlurredImage->Update();
		itkSobelImage = gaussian.EdgeDetection(itkBlurredImage);
		itkSobelImage->Update();



		adapter.ITKtoSTK(itkSobelImage, myResultResize2);

		stk::IO<float> imageIO;
		//imageIO.WriteImage( myImage, output );



		stk::ImageHistogram<TH2F, float> myImageHistogram;

		myImageHistogram.SetYAxisTitle("ADC");
		myImageHistogram.SetYAxisLog(false);
		myImageHistogram.SetNumberOfYBins(1024);
		myImageHistogram.SetYLimits( -0.5, 1023.5 );
		myImageHistogram.SetXAxisTitle( "Row" );
		myImageHistogram.SetNumberOfXBins( 1024);
		myImageHistogram.SetXLimits(  -0.5, 1023.5  );
		myImageHistogram.SetTitle("4x4 MLC Edge Detection");

		myImageHistogram.SetGridY(false);
		myImageHistogram.SetGridX(false);

		myImageHistogram.SetStatBoxOptions(0);
		myImageHistogram.SetOutputFileNameAndPath("blurred");
		myImageHistogram.SetOutputFileType( stk::FileType::PNG );

		myImageHistogram.Generate2DHistogram(myResultResize2, "Sobel");
		//myImageHistogram.GenerateXSlice(512);
		myImageHistogram.ExportRoot();



		myImageHistogram.SaveHistogram();



	return 1;
}
