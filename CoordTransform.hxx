// Copyright (C) 2003 Vivien Mallet
//
// This file is part of AtmoData library.
// AtmoData library is a tool for data processing in atmospheric
// sciences.
// 
// AtmoData is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// 
// AtmoData is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License (file "license") for more details.
//
// For more information, please see the AtmoData home page:
//     http://spacetown.free.fr/lib/atmodata


#ifndef ATMODATA_FILE_COORDTRANSFORM_HXX

#include <cmath>

namespace AtmoData
{
  
  //! Coordinate transformation from Lambert azimuthal equal area
  //! to longitude/latitude.
  template<class T>
  class LaeaToLonlat
  {
  protected:
    //! Earth radius.
    const T Earth_radius_;
    //! pi.
    const T pi_;
    //! Numerical limit assumed to be 0.
    const T limit_;

    //! Latitude of origin.
    T lat_origin_;
    //! Longitude of origin.
    T lon_origin_;

  public:
    LaeaToLonlat(T lon_origin, T lat_origin)  throw();
    void operator() (const T x, const T y,
		     T& lon, T& lat);
  };


}  // namespace AtmoData.

#define ATMODATA_FILE_COORDTRANSFORM_HXX
#endif