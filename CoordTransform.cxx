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


#ifndef ATMODATA_FILE_COORDTRANSFORM_CXX

#include "CoordTransform.hxx"

namespace AtmoData
{


  //! Default constructor.
  template <class T>
  LaeaToLonlat<T>::LaeaToLonlat(T lon_origin, T lat_origin)  throw():
    pi_(3.14159265358979323846264),
    lon_origin_(lon_origin / 180. * pi_),
    lat_origin_(lat_origin / 180. * pi_),
    Earth_radius_(6370997.), limit_(1.e-15)
  {

  }


  //! Convertion operator.
  /*!
    \param x_ abscissa in Lambert azimuthal equal area.
    \param y_ ordinate in Lambert azimuthal equal area.
    \param lon longitude (output).
    \param lat latitude (input).
  */
  template <class T>
  void LaeaToLonlat<T>::operator() (const T x_, const T y_,
				    T& lon, T& lat)
  {

    T rho, z, temp, cos_z, sin_z;
    T x(x_/Earth_radius_), y(y_/Earth_radius_);

    rho = sqrt(x * x + y * y);
    z = 2.0 * asin(rho / 2.0);
    cos_z = cos(z);
    sin_z = sin(z);

    if ((rho = hypot(x, y)) < limit_)
      {
	lon = lon_origin_;
	lat = lat_origin_;
	return;
      }

    x *= sin_z;
    T ab = cos_z * sin(lat_origin_) + y * sin_z * cos(lat_origin_) / rho;
    y = rho * cos(lat_origin_) * cos_z - y * sin(lat_origin_) * sin_z;

    lon = atan2(x, y);
    lat = asin(ab);
    lon += lon_origin_;
    if (fabs(lon) > pi_)
      {
	lon += pi_;
	lon -= 2.0 * pi_ * floor(lon / (2.0*pi_));
	lon -= pi_;
      }

    lat = lat/pi_*180.;
    lon = lon/pi_*180.;

  }


}  // namespace AtmoData.

#define ATMODATA_FILE_COORDTRANSFORM_CXX
#endif
