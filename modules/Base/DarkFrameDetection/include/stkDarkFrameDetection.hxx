/*
 * stkDarkFrameDetection.hxx
 *
 *  Created on: 4 Aug 2015
 *      Author: lewish
 */

#ifndef __stkDarkFrameDetection_hxx
#define __stkDarkFrameDetection_hxx


#include <cmath>
#include "stkDarkFrameDetection.h"
#include "stkImageMinus.h"


namespace stk{

DarkFrameDetection::DarkFrameDetection(){}
DarkFrameDetection::~DarkFrameDetection(){}






template< typename T_PixelInputType >
void DarkFrameDetection::CheckDarkFramesBefore(std::shared_ptr  < stk::ImageStack<T_PixelInputType> > myImageStack, int &iFrames, float &darkFrames, int &darkFlag, const float &pedAvg, const int &framesize){

	do {
		float tempPix=0;
		for(int iElements=0; iElements<framesize;iElements++)
		{
			tempPix+=myImageStack->GetPixelAt(iFrames*framesize+iElements);
		}
		tempPix=tempPix/framesize;

		if	(tempPix>(pedAvg+200))
		{
			darkFrames=iFrames;
			darkFlag=1;
		}
		iFrames++;

	}while(darkFrames==0&&iFrames<myImageStack->NumberOfImageInStack()&&darkFlag!=1);
}


template< typename T_PixelInputType >
void DarkFrameDetection::CheckDarkFramesAfter(std::shared_ptr  < stk::ImageStack<T_PixelInputType> > myImageStack, int &iFrames, float &darkFramesAfter, int &darkFlag, const float &pedAvg, const int &framesize){

	do {
		float tempPix=0;
		for(int iElements=0; iElements<framesize;iElements++)
		{
			tempPix+=myImageStack->GetPixelAt(iFrames*framesize+iElements);
		}
		tempPix=tempPix/framesize;

		if	(tempPix<(pedAvg+200))
		{
			darkFramesAfter=iFrames;

		}

		iFrames++;

	}while(darkFramesAfter==0&&iFrames<myImageStack->NumberOfImageInStack());
}





template< typename T_PixelOutputType >
void DarkFrameDetection::DetectDarkFrames( std::shared_ptr < stk::Image<T_PixelOutputType> > result, const std::shared_ptr < stk::Image<T_PixelOutputType> > pedestalImage, const float &pedAvg, const std::string &filePath, const std::string &fileNameAndFormat, const int &rows, const int &cols, const int  &numOfFrames, const int  &startingframe){
	stk::ImageMinus myMinus;
	int NumOfImageStacks(0);
	int framesize=rows*cols;
	double loops(0), dec(0), newframenumber(0);
	newframenumber=numOfFrames;
	dec = modf(newframenumber/100, &loops);


	std::cout<<"loops: "<<loops<<" "<<newframenumber<<std::endl;
	if( dec == 0 ){
		NumOfImageStacks = loops;
	}
	else{
		NumOfImageStacks = loops +1;
	}

	std::shared_ptr< stk::ImageStack<unsigned short> > myImageStack( new stk::ImageStack<unsigned short> );
	stk::IOImageStack<unsigned short> myIO;


	int darkFlag(0);
	for ( int iter = 0; iter<NumOfImageStacks;iter++){
		std::cout<<"iter: "<<iter<<std::endl;
		if (iter==loops)
		{

			std::cout<<"Iter is equal to loops"<<std::endl;
			if(NumOfImageStacks==1){
				myImageStack->Initialise( myIO.ReadImageStack( filePath, fileNameAndFormat, startingframe, (100*dec), framesize ), 100*dec	, rows, cols );
			}
			else
			{
				myImageStack->Initialise( myIO.ReadImageStack( filePath, fileNameAndFormat, iter*100, (100*dec), framesize ), 100*dec	, rows, cols );
			}
			//Calculating the number of dark frames in front of the data.

			float darkFrames=0;
			float darkFramesAfter=0;


			int iFrames=0;
			CheckDarkFramesBefore(myImageStack, iFrames, darkFrames, darkFlag, pedAvg, framesize);
			std::cout<<" Before dark done"<<std::endl;

			//Calculating the number of dark frames after the data.
			iFrames=darkFrames;
			CheckDarkFramesAfter(myImageStack, iFrames, darkFramesAfter, darkFlag, pedAvg, framesize);
			if(darkFramesAfter==0){
				darkFramesAfter=dec*100;
			}
			//Remove pedestal from raw stack, sum stack into image.
			myMinus.MinusImage(myImageStack, pedestalImage,result, darkFrames, darkFramesAfter);



		}
		else
		{

			myImageStack->Initialise( myIO.ReadImageStack( filePath, fileNameAndFormat, iter*100, (100), framesize ), 100, rows, cols );

			float darkFrames=0;
			float darkFramesAfter=0;
			int iFrames=0;

			CheckDarkFramesBefore(myImageStack, iFrames, darkFrames, darkFlag, pedAvg, framesize);
			//Calculating the number of dark frames after the data.
			iFrames=darkFrames;
			CheckDarkFramesAfter(myImageStack, iFrames, darkFramesAfter, darkFlag, pedAvg, framesize);
			if(darkFramesAfter==0){
				darkFramesAfter=100;
			}

			//Remove pedestal from raw stack, sum stack into image.
			myMinus.MinusImage(myImageStack, pedestalImage, result, darkFrames, darkFramesAfter);
			//myDivider.DivideImage(myResult, static_cast<float>(darkFramesAfter-darkFrames) );
		}

	}
}



template< typename T_PixelOutputType >
void DarkFrameDetection::DetectDarkFrames( std::shared_ptr < stk::Image<T_PixelOutputType> > result, const float &pedAvg, const std::string &filePath, const std::string &fileNameAndFormat, const int &rows, const int &cols, const int  &numOfFrames, const int  &startingframe, int &framesUsed ){


	stk::ImageSum mySummer;
	int NumOfImageStacks(0);
	int framesize=rows*cols;
	double loops(0), dec(0), newframenumber(0);
	newframenumber=numOfFrames;
	dec = modf(newframenumber/100, &loops);
	std::cout<<"loops: "<<loops<<" "<<newframenumber<<std::endl;
	if( dec == 0 ){
		NumOfImageStacks = loops;
	}
	else{
		NumOfImageStacks = loops +1;
	}

	std::shared_ptr< stk::ImageStack<unsigned short> > myImageStack( new stk::ImageStack<unsigned short> );
	stk::IOImageStack<unsigned short> myIO;


	int darkFlag(0);
	for ( int iter = 0; iter<NumOfImageStacks;iter++){
		std::cout<<"iter: "<<iter<<std::endl;
		if (iter==loops)
		{

			std::cout<<"Iter is equal to loops"<<std::endl;
			if(NumOfImageStacks==1){
				myImageStack->Initialise( myIO.ReadImageStack( filePath, fileNameAndFormat, startingframe, (100*dec), framesize ), 100*dec	, rows, cols );
			}
			else
			{
				myImageStack->Initialise( myIO.ReadImageStack( filePath, fileNameAndFormat, iter*100, (100*dec), framesize ), 100*dec	, rows, cols );
			}
			//Calculating the number of dark frames in front of the data.

			float darkFrames=0;
			float darkFramesAfter=0;


			std::cout<<"Starting first loop"<<std::endl;
			std::cout<<myImageStack->NumberOfImageInStack()<<std::endl;
			int iFrames=0;
			do {
				float tempPix=0;
				for(int iElements=0; iElements<framesize;iElements++)
				{
					tempPix+=myImageStack->GetPixelAt(iFrames*framesize+iElements);
				}
				tempPix=tempPix/framesize;
				std::cout<<"done first bit"<<std::endl;
				if	(tempPix>(pedAvg+200))
				{

					darkFrames=iFrames;
				}
				iFrames++;
				std::cout<<darkFrames<<" "<<darkFlag<<std::endl;
			}while(darkFrames==0&&iFrames<myImageStack->NumberOfImageInStack()&&darkFlag!=1);

			std::cout<<" Before dark done"<<std::endl;

			//Calculating the number of dark frames after the data.
			iFrames=darkFrames;
			do {
				float tempPix=0;
				for(int iElements=0; iElements<framesize;iElements++)
				{
					tempPix+=myImageStack->GetPixelAt(iFrames*framesize+iElements);
				}
				tempPix=tempPix/framesize;

				if	(tempPix<(pedAvg+200))
				{
					darkFramesAfter=iFrames;

				}

				iFrames++;

			}while(darkFramesAfter==0&&iFrames<myImageStack->NumberOfImageInStack());
			if(darkFramesAfter==0){
				darkFramesAfter=dec*100;
			}
			std::cout<<"Dark frames excluded."<<std::endl;
			std::cout<<darkFrames<<" "<<darkFramesAfter<<std::endl;
			//Remove pedestal from raw stack, sum stack into image.
			//
			mySummer.SumImageStack(myImageStack, result, darkFrames, darkFramesAfter);
			framesUsed+=(darkFramesAfter-darkFrames);
			//myDivider.DivideImage(myResult, static_cast<float>(darkFramesAfter-darkFrames) );
			//

		}
		else
		{


			std::cout<<iter*100<<" "<<framesize<<" "<<rows<<cols<<std::endl;

			myImageStack->Initialise( myIO.ReadImageStack( filePath, fileNameAndFormat, iter*100, (100), framesize ), 100, rows, cols );


			float darkFrames=0;
			float darkFramesAfter=0;


			int iFrames=0;
			do {
				float tempPix=0;
				for(int iElements=0; iElements<framesize;iElements++)
				{
					tempPix+=myImageStack->GetPixelAt(iFrames*framesize+iElements);
				}
				tempPix=tempPix/framesize;

				if	(tempPix>(pedAvg+200))
				{

					darkFrames=iFrames;
					darkFlag=1;
				}
				iFrames++;

			}while(darkFrames==0&&iFrames<myImageStack->NumberOfImageInStack()&&darkFlag!=1);



			//Calculating the number of dark frames after the data.
			iFrames=darkFrames;
			do {
				float tempPix=0;
				for(int iElements=0; iElements<framesize;iElements++)
				{
					tempPix+=myImageStack->GetPixelAt(iFrames*framesize+iElements);
				}
				tempPix=tempPix/framesize;

				if	(tempPix<(pedAvg+200))
				{
					darkFramesAfter=iFrames;
				}

				iFrames++;


			}while(darkFramesAfter==0&&iFrames<myImageStack->NumberOfImageInStack());
			if(darkFramesAfter==0){
				darkFramesAfter=100;
			}
			std::cout<<darkFrames<<" "<<darkFramesAfter<<std::endl;
			std::cout<<"Dark frames excluded."<<std::endl;
			//Remove pedestal from raw stack, sum stack into image.
			mySummer.SumImageStack(myImageStack, result, darkFrames, darkFramesAfter);
			framesUsed+=(darkFramesAfter-darkFrames);
			//myDivider.DivideImage(myResult, static_cast<float>(darkFramesAfter-darkFrames) );
		}

	}
}

}


#endif /*__stkDarkFrameDetection_h */
