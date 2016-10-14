/*============================
 * stkImageSum.hxx
 *
 *  Created on: 8 Jun 2015
 *      Author: phrfp
 *=============================*/

#ifndef __stkImageSum_hxx
#define __stkImageSum_hxx

#include <algorithm>
#include <functional>
#include "stkImageSum.h"

namespace stk{


ImageSum::ImageSum(){}


ImageSum::~ImageSum(){}



template< typename T_PixelInputType, typename T_PixelOutputType>
void ImageSum::SumImageStack( const std::shared_ptr< stk::ImageStack<T_PixelInputType> > imageStack, std::shared_ptr < stk::Image<T_PixelOutputType> > result ){

	typename std::vector<T_PixelInputType>::iterator itFrameStart;
	typename std::vector<T_PixelInputType>::iterator itFrameEnd;
	typename std::vector<T_PixelOutputType>::iterator itResult;

	for(int iFrames=0; iFrames < imageStack->NumberOfImageInStack(); iFrames++)

	{
		itFrameStart = imageStack->StartOfFrame(iFrames);
		itFrameEnd = imageStack->EndOfFrame(iFrames);
		itResult = result->StartImage();
		std::transform( itFrameStart, itFrameEnd, itResult, itResult, std::plus<T_PixelOutputType>() );
		//std::cout<<iFrames<<std::endl;
	}
}

template< typename T_PixelInputType, typename T_PixelOutputType>
	void ImageSum::SumImageStack( const std::shared_ptr< stk::ImageStack<T_PixelInputType> > imageStack, std::shared_ptr < stk::Image<T_PixelOutputType> > result, const int &darkFrames, const int &darkFramesAfter ){

	for(int iFrames=darkFrames; iFrames<darkFramesAfter; iFrames++)
			{

		for(int iElements=0; iElements<result->NumberOfPixels();iElements++){
			T_PixelOutputType temp(0);
			temp = imageStack->GetPixelAt(iFrames*result->NumberOfPixels()+iElements)+result->GetPixelAt(iElements);


			result->SetPixelAt(iElements,temp );
		}

			}
}


}
#endif /* __stkImageSum_hxx */
