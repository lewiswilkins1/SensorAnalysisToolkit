/*
 * stkDarkFrameDetection.h
 *
 *  Created on: 4 Aug 2015
 *      Author: lewish
 */

#ifndef __stkDarkFrameDetection_h
#define __stkDarkFrameDetection_h

#include <memory>
#include <iterator>
#include <vector>

#include "stkImage.h"
#include "stkImageStack.h"



namespace stk{


class DarkFrameDetection{


public:
	DarkFrameDetection();
	virtual ~DarkFrameDetection();

	template< typename T_PixelInputType >
	void CheckDarkFramesBefore(std::shared_ptr  < stk::ImageStack<T_PixelInputType> > myImageStack, int &iFrames, float &darkFrames, int &darkFlag, const float &pedAvg,const int &framesize);

	template< typename T_PixelInputType >
	void CheckDarkFramesAfter(std::shared_ptr  < stk::ImageStack<T_PixelInputType> > myImageStack, int &iFrames, float &darkFramesAfter, int &darkFlag, const float &pedAvg,const int &framesize);

	template< typename T_PixelOutputType >
	void DetectDarkFrames( std::shared_ptr < stk::Image<T_PixelOutputType> > result, const std::shared_ptr < stk::Image<T_PixelOutputType> > pedestalImage, const float &pedAvg, const std::string &filePath, const std::string &fileNameAndFormat, const int &rows, const int &cols, const int  &numOfFrames, const int  &startingframe);

	template< typename T_PixelOutputType >
	void DetectDarkFrames( std::shared_ptr < stk::Image<T_PixelOutputType> > result, const float &pedAvg, const std::string &filePath, const std::string &fileNameAndFormat, const int &rows, const int &cols, const int  &numOfFrames, const int  &startingframe, int &framesUsed  );

};


}

#include "stkDarkFrameDetection.hxx"


#endif /* __stkDarkFrameDetection_h */
