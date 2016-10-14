/*
 * SliceComparison.cxx
 *
 *  Created on: 31 Jul 2015
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
//#include "stkImageFlatFieldCorrection.h"
#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "stkITKTools.h"
#include "itkImageFileWriter.h"
#include "stkITKAdapter.h"
#include "stkITKTools.h"

int main(int argc, const char** argv){

	std::string filePathandfilename, otherfile;
	int startingframe, numOfFrames, rows, cols, framesize;

	std::stringstream inputs;
	for(int iArg=1; iArg < argc; iArg++)
	{
		inputs << argv[iArg] << ' ';//get the arguments
	}

	inputs >> filePathandfilename >> otherfile >> rows >> cols ;//write the inputs
	framesize = cols*rows;


	stk::IO<float> myIO;
	std::shared_ptr< stk::Image<float> > myImage( new stk::Image<float> );
	myImage->Initialise( myIO.ReadImage( filePathandfilename ), rows, cols );
	stk::IO<float> myIO2;
	std::shared_ptr< stk::Image<float> > myImage2( new stk::Image<float> );
	myImage2->Initialise( myIO.ReadImage( otherfile ), rows, cols );

	//Image to hold resized result
	std::shared_ptr<stk::Image<float> > myResultResize (new stk::Image<float>(1024,1024));
	std::shared_ptr<stk::Image<float> > myResultResize2 (new stk::Image<float>(1024,1024));

	stk::ImageResize mySize;

	//Image resizing
	mySize.ResizeImage(myImage , myResultResize);
	mySize.ResizeImage(myImage2 , myResultResize2);

	typename std::vector<float>::iterator imageStart;
	typename std::vector<float>::iterator imageEnd;
	typename std::vector<float>::iterator imageStart2;
	typename std::vector<float>::iterator imageEnd2;
	imageStart=myResultResize->StartImage();
	imageEnd=myResultResize->EndImage();
	imageStart2=myResultResize2->StartImage();
	imageEnd2=myResultResize2->EndImage();

	float min(0), max(0), min2(0), max2(0);
	min = *std::min_element(imageStart,imageEnd);
	max = *std::max_element(imageStart,imageEnd);

	min2 = *std::min_element(imageStart2,imageEnd2);
	max2 = *std::max_element(imageStart2,imageEnd2);
	std::cout<<"Max: "<<max<<" Min: "<<min<<std::endl;

	for(int iElements=0; iElements<1024*1024; iElements++){
		float Temp, temp2;
		Temp = (myResultResize->GetPixelAt(iElements))/(max);
		temp2 = (myResultResize2->GetPixelAt(iElements))/(max2);
		myResultResize->SetPixelAt(iElements, Temp);
		myResultResize2->SetPixelAt(iElements, temp2);

	}






	typedef itk::Image< float, 2 >  ImageType;
	ImageType::Pointer itkImage = ImageType::New();
	ImageType::Pointer itkBlurredImage = ImageType::New();
	ImageType::Pointer itkImage2 = ImageType::New();
	ImageType::Pointer itkBlurredImage2 = ImageType::New();
	ImageType::Pointer itkTranslated = ImageType::New();



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
	itkBlurredImage2->SetRegions(region);
	itkBlurredImage2->Allocate();
	itkTranslated->SetRegions(region);
	itkTranslated->Allocate();
	stk::ITKAdapter adapter;
	adapter.STKtoITK(itkImage, myResultResize, 1024, 1024);
	adapter.STKtoITK(itkImage2, myResultResize2, 1024, 1024);
	stk::ITKTools gaussian;

	itkBlurredImage = gaussian.GaussianBlur(itkImage);
	itkBlurredImage->Update();

	itkBlurredImage2 = gaussian.GaussianBlur(itkImage2);
	itkBlurredImage2->Update();


	itkTranslated = gaussian.ImageRegistration(itkBlurredImage2, itkBlurredImage);

	itkTranslated->Update();



	adapter.ITKtoSTK(itkTranslated, myResultResize2);
	adapter.ITKtoSTK(itkBlurredImage2, myResultResize);




	stk::ImageHistogram<TH2F, float> myImageHistogram;

	myImageHistogram.SetYAxisTitle("Normalised Counts");
	myImageHistogram.SetYAxisLog(false);
	myImageHistogram.SetNumberOfYBins(1024);
	myImageHistogram.SetYLimits( -0.5, 1023.5 );
	myImageHistogram.SetXAxisTitle( "Col" );
	myImageHistogram.SetNumberOfXBins( 1024);
	myImageHistogram.SetXLimits(  -0.5, 1023.5  );

	myImageHistogram.SetGridY(false);
	myImageHistogram.SetGridX(false);

	myImageHistogram.SetStatBoxOptions(0);
	myImageHistogram.SetOutputFileNameAndPath("5x5_BU_vs_NoBU_Y_Slice");
	myImageHistogram.SetOutputFileType( stk::FileType::PNG );

	myImageHistogram.Generate2DHistogram( myResultResize2);
	myImageHistogram.Generate2DHistogram( myResultResize);
	//myImageHistogram.GenerateYSlice(512);
	myImageHistogram.SetTitle("5x5 BU vs NoBU Y Slice");
	myImageHistogram.SaveHistogram();


	stk::ImageHistogram<TH2F, float> myImageHistogram2;

	myImageHistogram2.SetYAxisTitle("Normalised Counts");
	myImageHistogram2.SetYAxisLog(false);
	myImageHistogram2.SetNumberOfYBins(1024);
	myImageHistogram2.SetYLimits( -0.5, 1023.5 );
	myImageHistogram2.SetXAxisTitle( "Col" );
	myImageHistogram2.SetNumberOfXBins( 1024);
	myImageHistogram2.SetXLimits(  -0.5, 1023.5  );

	myImageHistogram2.SetGridY(false);
	myImageHistogram2.SetGridX(false);

	myImageHistogram2.SetStatBoxOptions(0);
	myImageHistogram2.SetOutputFileNameAndPath("5x5_BU_vs_NoBU_X_Slice");
	myImageHistogram2.SetOutputFileType( stk::FileType::PNG );

	myImageHistogram2.Generate2DHistogram( myResultResize2);
	myImageHistogram2.Generate2DHistogram( myResultResize);

	//myImageHistogram2.GenerateXSlice(512);
	myImageHistogram2.SetTitle("5x5 BU vs NoBU X Slice");
	myImageHistogram2.SaveHistogram();

}





