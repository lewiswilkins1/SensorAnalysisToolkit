/*
 * stkITKAdapter.h
 *
 *  Created on: 4 Aug 2015
 *      Author: lewish
 */

#ifndef __stkITKAdapter_h
#define __stkITKAdapter_h


#include <memory>
#include <iterator>
#include <vector>

#include "stkImage.h"
#include "stkImageStack.h"
#include "itkImage.h"
#include "itkImageRegionIterator.h"
#include "itkImageFileWriter.h"
//#include "itkImageBase.h"


namespace stk{



class ITKAdapter{

public:

	ITKAdapter();
	virtual ~ITKAdapter();

	template<typename T_Pixeltype>
	void STKtoITK(itk::Image< float, 2 >::Pointer itkOutput, const std::shared_ptr< stk::Image<T_Pixeltype> > inputImage, const int &rows, const int &cols);


	template<typename T_Pixeltype>
	void ITKtoSTK(itk::Image< float, 2 >::Pointer itkInput, const std::shared_ptr< stk::Image<T_Pixeltype> > outputImage);



};

}

#include "stkITKAdapter.hxx"

#endif /*__stkImage_h */
