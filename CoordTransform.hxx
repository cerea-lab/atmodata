// Copyright (C) 2003-2004 CEREA
//     Author: Vivien Mallet
//
// CEREA (http://www.enpc.fr/cerea) is a joint laboratory of
// ENPC (http://www.enpc.fr) and EDF R&D (http://www.edf.fr).
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


  //! Coordinate transformation from indices of the MM5 grid
  //! in Lambert conic conformal to longitude/latitude.
  template<class T>
  class MM5LccIndToLonlat
  {
  protected:
    //! Coarse domain grid dimension in East-West direction.
    const int jmx_;
    //! Coarse domain grid dimension in North-South direction.
    const int imx_;
    //! East-West location in the coarse domain of the South-West corner.
    const double jx_;
    //! North-South location in the coarse domain of the South-West corner.
    const double ix_;
    //! Coarse domain center latitude (degree).
    const double phic_;
    //! Coarse domain center longitude (degree).
    const double lambdac_;
    //! True latitude #1 (degree).
    const double phi1_;
    //! True latitude #2 (degree).
    const double phi2_;
    //! Grid distance (meters) of the current domain.
    const double ds_;
    //! Domain grid size ratio with respect to coarse domain.
    const int ratio_;
    //! Earth radius.
    const double Earth_radius_;
    //! pi.
    const double pi_;

  public:
    MM5LccIndToLonlat(int jmx, int imx, double jx, double ix, double phic, double lambdac,
		      double phi1, double phi2, double ds, int ratio)  throw();
    void operator() (const T j, const T i,
		     T& lon, T& lat);
  };


  //! Coordinate transformation from longitude/latitude to
  //! indices of the MM5 grid in Lambert conic conformal.
  template<class T>
  class LonlatToMM5LccInd
  {
  protected:
    //! Coarse domain grid dimension in East-West direction.
    const int jmx_;
    //! Coarse domain grid dimension in North-South direction.
    const int imx_;
    //! East-West location in the coarse domain of the South-West corner.
    const double jx_;
    //! North-South location in the coarse domain of the South-West corner.
    const double ix_;
    //! Coarse domain center latitude (degree).
    const double phic_;
    //! Coarse domain center longitude (degree).
    const double lambdac_;
    //! True latitude #1 (degree).
    const double phi1_;
    //! True latitude #2 (degree).
    const double phi2_;
    //! Grid distance (meters) of the coarse domain.
    const double ds0_;
    //! Domain grid size ratio with respect to coarse domain.
    const int ratio_;
    //! Earth radius.
    const double Earth_radius_;
    //! pi.
    const double pi_;

  public:
    LonlatToMM5LccInd(int jmx, int imx, double jx, double ix, double phic, double lambdac,
		      double phi1, double phi2, double ds0, int ratio)  throw();
    void operator() (const T lon, const T lat,
		     T& j, T& i);
  };


  //! Coordinate transformation from indices of the MM5 grid
  //! in Mercator coordinates to longitude/latitude.
  template<class T>
  class MM5MercIndToLonlat
  {
  protected:
    //! Coarse domain grid dimension in East-West direction.
    const int jmx_;
    //! Coarse domain grid dimension in North-South direction.
    const int imx_;
    //! East-West location in the coarse domain of the South-West corner.
    const double jx_;
    //! North-South location in the coarse domain of the South-West corner.
    const double ix_;
    //! Coarse domain center latitude (degree).
    const double phic_;
    //! Coarse domain center longitude (degree).
    const double lambdac_;
    //! True latitude #1 (degree).
    const double phi1_;
    //! Grid distance (meters) of the current domain.
    const double ds_;
    //! Domain grid size ratio with respect to coarse domain.
    const int ratio_;
    //! Earth radius.
    const double Earth_radius_;
    //! pi.
    const double pi_;

  public:
    MM5MercIndToLonlat(int jmx, int imx, double jx, double ix, double phic, double lambdac,
		       double phi1, double ds, int ratio)  throw();
    void operator() (const T j, const T i,
		     T& lon, T& lat);
  };


  //! Coordinate transformation from longitude/latitude to
  //! indices of the MM5 grid in Mercator coordinates.
  template<class T>
  class LonlatToMM5MercInd
  {
  protected:
    //! Coarse domain grid dimension in East-West direction.
    const int jmx_;
    //! Coarse domain grid dimension in North-South direction.
    const int imx_;
    //! East-West location in the coarse domain of the South-West corner.
    const double jx_;
    //! North-South location in the coarse domain of the South-West corner.
    const double ix_;
    //! Coarse domain center latitude (degree).
    const double phic_;
    //! Coarse domain center longitude (degree).
    const double lambdac_;
    //! True latitude #1 (degree).
    const double phi1_;
    //! Grid distance (meters) of the coarse domain.
    const double ds0_;
    //! Domain grid size ratio with respect to coarse domain.
    const int ratio_;
    //! Earth radius.
    const double Earth_radius_;
    //! pi.
    const double pi_;

  public:
    LonlatToMM5MercInd(int jmx, int imx, double jx, double ix, double phic, double lambdac,
		       double phi1, double ds0, int ratio)  throw();
    void operator() (const T lon, const T lat,
		     T& j, T& i);
  };


  //! Coordinate transformation from indices of the MM5 grid
  //! in polar stereographic coordinates to longitude/latitude.
  template<class T>
  class MM5StereIndToLonlat
  {
  protected:
    //! Coarse domain grid dimension in East-West direction.
    const int jmx_;
    //! Coarse domain grid dimension in North-South direction.
    const int imx_;
    //! East-West location in the coarse domain of the South-West corner.
    const double jx_;
    //! North-South location in the coarse domain of the South-West corner.
    const double ix_;
    //! Coarse domain center latitude (degree).
    const double phic_;
    //! Coarse domain center longitude (degree).
    const double lambdac_;
    //! True latitude #1 (degree).
    const double phi1_;
    //! Grid distance (meters) of the current domain.
    const double ds_;
    //! Domain grid size ratio with respect to coarse domain.
    const int ratio_;
    //! Earth radius.
    const double Earth_radius_;
    //! pi.
    const double pi_;

  public:
    MM5StereIndToLonlat(int jmx, int imx, double jx, double ix, double phic, double lambdac,
			double phi1, double ds, int ratio)  throw();
    void operator() (const T j, const T i,
		     T& lon, T& lat);
  };


  //! Coordinate transformation from longitude/latitude to
  //! indices of the MM5 grid in polar stereographic coordinates.
  template<class T>
  class LonlatToMM5StereInd
  {
  protected:
    //! Coarse domain grid dimension in East-West direction.
    const int jmx_;
    //! Coarse domain grid dimension in North-South direction.
    const int imx_;
    //! East-West location in the coarse domain of the South-West corner.
    const double jx_;
    //! North-South location in the coarse domain of the South-West corner.
    const double ix_;
    //! Coarse domain center latitude (degree).
    const double phic_;
    //! Coarse domain center longitude (degree).
    const double lambdac_;
    //! True latitude #1 (degree).
    const double phi1_;
    //! Grid distance (meters) of the coarse domain.
    const double ds0_;
    //! Domain grid size ratio with respect to coarse domain.
    const int ratio_;
    //! Earth radius.
    const double Earth_radius_;
    //! pi.
    const double pi_;

  public:
    LonlatToMM5StereInd(int jmx, int imx, double jx, double ix, double phic, double lambdac,
			double phi1, double ds0, int ratio)  throw();
    void operator() (const T lon, const T lat,
		     T& j, T& i);
  };


}  // namespace AtmoData.

#define ATMODATA_FILE_COORDTRANSFORM_HXX
#endif
