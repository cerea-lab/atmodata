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


#ifndef ATMODATA_FILE_PHOTOLYSIS_CXX

#include "Photolysis.hxx"

namespace AtmoData
{


  // The following function is derived from the Fortran (77) subroutine
  // zenith found in TUV (Tropospheric Ultraviolet & Visible radiation model).
  // This subroutine was provided under the GNU General Public License
  // under the following copyright:
  // Copyright (C) 1994,95,96  University Corporation for Atmospheric Research.
  /*! Calculates solar zenith angle for a given time and location.
    Calculation is based on equations given in:  Paltridge and Platt,
    Radiative Processes in Meteorology and Climatology, Elsevier, pp. 62,63, 1976.
    Fourier coefficients originally from:  Spencer, J.W., 1971, Fourier
    series representation of the position of the sun, Search, 2:172.
    Note:  This approximate program does not account for changes from year
    to year.
    \param lat latitude of location (degrees).
    \param long longitude of location (degrees).
    \param idate date in the form YYMMDD.
    \param ut local time in decimal UT (e.g., 16.25 means 15 minutes after 4 pm).
    \return solar zenith angle (degrees).
  */
  template<class T>
  T Zenith(T lat, T lon, int idate, T ut)
  {

    T azim, zen;

    T lbut,lzut;
    T rlt;
    T d, tz, rdecl, eqr, eqh, zpt;
    T csz, zr, caz, raz ;
    T sintz, costz, sin2tz, cos2tz, sin3tz, cos3tz;

    int iiyear, imth, iday, ijd;
    int imn[] = {31, 28, 31, 30, 31, 30,
		 31, 31, 30, 31, 30, 31};
  
    int i;

    const T pi = 3.1415926535898;
    const T dr = double(pi) / double(180.);

    // Converts to radians.
    rlt = lat*dr;

    // Parse date.
    iiyear = idate / 10000;
    imth = (idate - iiyear * 10000) / 100;
    iday = idate - iiyear * 10000 - imth * 100;

    // Identifies and corrects leap years.
    if (iiyear%4 == 0)
      imn[1] = 29;
    else
      imn[1] = 28;

    // Computes current (Julian) day of year IJD = 1 to 365.
    ijd = 0;
    for (i=0; i<imth-1; i++)
      ijd += imn[i];
    ijd += iday;

    // Calculates decimal Julian day from start of year.
    d = T(ijd - 1) + ut / 24.;

    // Equation 3.8 for "day-angle".
    tz = 2.* pi * d / 365.;

    // Calculates sine and cosine from addition theoremes for 
    // better performance;  the computation of sin2tz,
    // sin3tz, cos2tz and cos3tz is about 5-6 times faster
    // than the evaluation of the intrinsic functions.
    //
    // It is sin(x+y) = sin(x)*cos(y)+cos(x)*sin(y)
    // and   cos(x+y) = cos(x)*cos(y)-sin(x)*sin(y)
    //
    // sintz  = sin(tz)      costz  = cos(tz)
    // sin2tz = sin(2.*tz)   cos2tz = sin(2.*tz)
    // sin3tz = sin(3.*tz)   cos3tz = cos(3.*tz)

    sintz = sin(tz);
    costz = cos(tz);
    sin2tz = 2. * sintz * costz;
    cos2tz = costz * costz - sintz * sintz;
    sin3tz = sintz * cos2tz + costz * sin2tz;
    cos3tz = costz * cos2tz - sintz * sin2tz;

    // Equation 3.7 for declination in radians.
    rdecl = 0.006918 - 0.399912 * costz  + 0.070257 * sintz
      - 0.006758 * cos2tz + 0.000907 * sin2tz
      - 0.002697 * cos3tz + 0.001480 * sin3tz;

    // Equation 3.11 for Equation of time in radians.
    eqr   = 0.000075 + 0.001868 * costz  - 0.032077 * sintz
      - 0.014615 * cos2tz - 0.040849 * sin2tz;

    // Converts equation of time to hours.
    eqh = eqr * 24. / (2. * pi);

    // Calculates local hour angle (hours).
    lbut = 12. - eqh - lon* 24. / 360.;

    // Converts to angle from UT.
    lzut = 15. * (ut - lbut);
    zpt = lzut * dr;

    // Equation 2.4 for cosine of zenith angle.
    csz = sin(rlt) * sin(rdecl) + cos(rlt) * cos(rdecl) * cos(zpt);
    zr = acos(csz);
    zen = zr / dr;

    // Calculates local solar azimuth.  
    // caz = (sin(rdecl) - sin(rlt) * cos(zr)) / (cos(rlt) * sin(zr));
    // raz = acos(caz);
    // azim = raz / dr;
  
    return zen;
  
  }


}  // namespace AtmoData.

#define ATMODATA_FILE_PHOTOLYSIS_CXX
#endif
