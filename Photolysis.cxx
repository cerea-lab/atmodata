// Copyright (C) 2003-2004 Vivien Mallet
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
    \param lon longitude of location (degrees).
    \param lat latitude of location (degrees).
    \param idate date in the form YYYYMMDD.
    \param ut local time in decimal UT (e.g., 16.25 means 15 minutes after 4 pm).
    \return solar zenith angle (degrees).
  */
  template<class T>
  T ZenithAngle(T lon, T lat, int idate, T ut)
  {

    T azim, zen;

    T lbut,lzut;
    T rlt;
    T d, tz, rdecl, eqr, eqh, zpt;
    T csz, zr, caz, raz ;
    T sintz, costz, sin2tz, cos2tz, sin3tz, cos3tz;

    int iiiiyear, imth, iday, ijd;
    int imn[] = {31, 28, 31, 30, 31, 30,
		 31, 31, 30, 31, 30, 31};
  
    int i;

    const T pi = 3.1415926535898;
    const T dr = double(pi) / double(180.);

    // Converts to radians.
    rlt = lat*dr;

    // Parse date.
    iiiiyear = idate / 10000;
    imth = (idate - iiiiyear * 10000) / 100;
    iday = idate - iiiiyear * 10000 - imth * 100;

    // Identifies and corrects leap years.
    if ( ( (iiiiyear % 4 == 0) && (iiiiyear % 100 != 0) )
	 || ( iiiiyear % 400 == 0) )
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


  //! Computes the cloud attenuation for photolysis rates.
  /*!
    \param Temperature temperature (K).
    \param Pressure pressure (Pa).
    \param Humidity specific humidity (unitless).
    \param LiquidWaterContent liquid water content (kg/m^3).
    \param MediumCloudiness medium cloudiness (in [0, 1]).
    \param HighCloudiness high cloudiness (in [0, 1]).
    \param CriticalRelativeHumidity function that returns the critical
    relative humidity as function of the altitude, the pressure and reference pressure.
    \param date date in the form YYYYMMDD.
    \param Attenuation (output) cloud attenuation coefficient.
  */
  template <class TT, class TP, class TH, class TL,
	    class TMC, class THC, class T, class TG>
  void Attenuation_LWC(Data<TT, 4, TG>& Temperature, Data<TP, 4, TG>& Pressure,
		       Data<TH, 4, TG>& Humidity, Data<TL, 4, TG>& LiquidWaterContent,
		       Data<TMC, 3, TG>& MediumCloudiness, Data<THC, 3, TG>& HighCloudiness,
		       T (CriticalRelativeHumidity)(const T&, const T&, const T&),
		       int date, Data<T, 4, TG>& Attenuation)
  {

    int h, k, j, i;
    int Nt(Attenuation.GetLength(0));
    int Nz(Attenuation.GetLength(1));
    int Ny(Attenuation.GetLength(2));
    int Nx(Attenuation.GetLength(3));
    

    // Index "0" and "1" refer to two contiguous levels.
    T rh0, rh1, rhc, dz, delta_z,
      lwc0, lwc1, lw, w, tau, opt_depth, s, tr;

    for (h=0; h<Nt; h++)
      for (j=0; j<Ny; j++)
	for (i=0; i<Nx; i++)
	  {

	    // Specific humidity to relative humidity.
	    s = 611. * pow(10., 7.5 * (Temperature(h, Nz-1, j, i) - 273.15)
			   / (Temperature(h, Nz-1, j, i) - 35.85));
	    rh0 = Pressure(h, Nz-1, j, i) * Humidity(h, Nz-1, j, i)
	      / (0.62197 + Humidity(h, Nz-1, j, i)) / s;
	    // kg/m^3 to g/m^3.
	    lwc0 = 1000. * LiquidWaterContent(h, Nz-1, j, i);

	    w = 0;
	    
	    for (k=Nz-1; k>=0; k--)
	      {

		if (k==Nz-1)
		  dz = Attenuation[1].Value(h, Nz-1, j, i)
		    - Attenuation[1].Value(h, Nz-2, j, i);
		else
		  dz = Attenuation[1].Value(h, k+1, j, i)
		    - Attenuation[1].Value(h, k, j, i);

		// Specific humidity to relative humidity.
		s = 611. * pow(10., 7.5 * (Temperature(h, k, j, i) - 273.15)
			       / (Temperature(h, k, j, i) - 35.85));
		rh1 = Pressure(h, k, j, i) * Humidity(h, k, j, i)
		  / (0.62197 + Humidity(h, k, j, i)) / s;
		// kg/m^3 to g/m^3.
		lwc1 = 1000. * LiquidWaterContent(h, k, j, i);

		// Critical relative humidity.
		rhc = CriticalRelativeHumidity(Attenuation[1].Value(h, k, j, i),
					       Pressure(h, k, j, i),
					       Pressure(h, 0, j, i));

		if ( (rh0>rhc) && (rh1>rhc) )  // In a cloud.
		  w += dz * (lwc0 + lwc1) / 2.0;
		else if ( (rh0>rhc) && (rh1<rhc) )  // Below a cloud.
		  {
		    delta_z = dz * (rh0 - rhc) / (rh0 - rh1);
		    w += lwc1 * delta_z + (lwc0 - lwc1) / dz * .5
		      * (2.0 * dz - delta_z) * delta_z;
		  }
		else if ( (rh0<rhc) && (rh1>rhc) )  // Above a cloud.
		  {
		    delta_z = dz * (rh1 - rhc) / (rh1 - rh0);
		    w += lwc1 * delta_z + (lwc0 - lwc1) / dz * .5
		      * delta_z * delta_z;
		  }
		// For the next level.
		rh0 = rh1;
		lwc0 = lwc1;

		// Computes liquid water path.
		if (w>0.)
		  lw = log10(w);
		else
		  lw = 0.;

		// Computes the cloud optical depth according to Stephens (1978).
		if (lw<=0.)
		  tau = 0.;
		else
		  tau = pow(10., 0.2633 + 1.7095 * log(lw));

		// Computes the cloud transmissivity.
		if (tau<5.)
		  tr = 1.;
		else
		  tr = (5. - exp(-tau)) / (4. + 0.42*tau);

		// Zenith angle.
		T cos_zenith_angle;
		cos_zenith_angle = cos( ZenithAngle(Attenuation[3].Value(h, k, j, i),
						    Attenuation[2].Value(h, k, j, i),
						    date, Attenuation[0].Value(h, k, j, i))
					* 0.0174532925199433 );
		cos_zenith_angle = abs(cos_zenith_angle);

		// Computes the attenuation coefficient.
		if (tr==1)
		  Attenuation(h, k, j, i) = 1.0
		    + (min(1.0, MediumCloudiness(h, j, i) + HighCloudiness(h, j, i)))
		    * (1.6 * cos_zenith_angle - 1.0);
		else
		  Attenuation(h, k, j, i) = 1.0
		    + (min(1.0, MediumCloudiness(h, j, i) + HighCloudiness(h, j, i)))
		    * ( (1.-tr) * cos_zenith_angle );

	      }

	  }
  }


}  // namespace AtmoData.

#define ATMODATA_FILE_PHOTOLYSIS_CXX
#endif
