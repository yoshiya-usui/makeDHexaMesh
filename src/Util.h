//--------------------------------------------------------------------------
// MIT License
//
// Copyright (c) 2021 Yoshiya Usui
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//--------------------------------------------------------------------------
#ifndef DBLDEF_UTIL
#define DBLDEF_UTIL

#include <iostream>
#include <string>
#include <algorithm>
#include <vector>

#include "CommonParameters.h"

namespace Util{

// Calculate length between inputted two coordinates
double calcLength( const CommonParameters::XY& coord0, const CommonParameters::XY& coord1 );

// Calculate length between inputted two coordinates
double calcLength( const CommonParameters::XYZ& coord0, const CommonParameters::XYZ& coord1 );

// Average two coordinates
CommonParameters::XYZ averageTwoCoords( const CommonParameters::XYZ& coord1, const CommonParameters::XYZ& coord2 );

// Check whether the location of the two points are same
bool isSameLocation( const CommonParameters::XYZ& coord1, const CommonParameters::XYZ& coord2 );

}

#endif
