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
    \param lat latitude (ouput).
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


  //! Default constructor.
  template <class T>
  MM5LccIndToLonlat<T>::MM5LccIndToLonlat(int jmx, int imx, double jx, double ix, double phic, double lambdac,
					  double phi1, double phi2, double ds, int ratio)  throw():
    jmx_(jmx), imx_(imx), jx_(jx), ix_(ix), phic_(phic), lambdac_(lambdac),
    phi1_(phi1), phi2_(phi2), ds_(ds), ratio_(ratio),
    pi_(3.14159265358979323846264),
    Earth_radius_(6370997.)
  {
    
  }


  //! Convertion operator.
  /*!
    \param j index of the MM5 grid along the East-West direction.
    \param i index of the MM5 grid along the North-South direction.
    \param lon longitude (output).
    \param lat latitude (output).
    \warning Indices order (in MM5) is confusing: j has to be provided first.
  */
  template <class T>
  void MM5LccIndToLonlat<T>::operator() (const T j, const T i,
					 T& lon, T& lat)
  {

    double ic0, jc0;
    double ic, jc;
    double conv;
    double aux1, aux2, kappa;
    double auxsig;
    double psi1;
    double auxc, yc;
    double x, y, R;
    double auxl;
    double lambdaprima;
    double aux, Rs;
    double auxij;
    double auxtan;

    ic0 = (imx_ + 1.0) / 2.0;
    jc0 = (jmx_ + 1.0) / 2.0;


    ic = (ic0 - ix_) * ratio_ + 0.5;
    jc = (jc0 - jx_) * ratio_ + 0.5;

    conv = 360.0 / (2.0 * pi_);

    aux1 = (45.0 - abs(phi1_) / 2.0) / conv;
    aux2 = (45.0 - abs(phi2_) / 2.0) / conv;
    kappa = ( log10(cos(phi1_/conv)) - log10(cos(phi2_/conv)) )
      / ( log10(tan(aux1)) - log10(tan(aux2)) );

    auxsig = phic_>0 ? 1.0 : -1.0;
    
    psi1 = auxsig * ( pi_ / 2.0 - abs(phi1_) / conv );

    auxc = (auxsig * 90.0 - phic_) / conv / 2.0;
    yc = - (Earth_radius_ / kappa) * sin(psi1)
      * pow(tan(auxc) / tan(psi1/2.0), kappa);

    x = (jc - j - 1) * ds_;
    y = (i + 1 - ic) * ds_ + yc;
    R = sqrt(x*x + y*y);

    auxl = tan(psi1/2.0)
	  * pow(auxsig * R * kappa / (Earth_radius_*sin(psi1)),
		  1.0/kappa);
    lat = auxsig * 90.0 - 2.0 * conv * atan(auxl);

    if (y == 0.0)
      {
	auxtan = x<0.0 ? (-pi_/2.0) : (pi_/2);
	lambdaprima = lambdac_ + conv / kappa * auxtan;
      }
    else
      lambdaprima = lambdac_ + conv / kappa * atan(x / (auxsig*y));

    if (lambdaprima < 180.0)
      lon = lambdaprima + 360.0;
    if ((-180.0 <= lambdaprima) && (lambdaprima <= 180.0))
      lon = lambdaprima;
    if (lambdaprima > 180.0)
      lon = lambdaprima - 360.0;

  }


  //! Default constructor.
  template <class T>
  LonlatToMM5LccInd<T>::LonlatToMM5LccInd(int jmx, int imx, double jx, double ix, double phic, double lambdac,
					  double phi1, double phi2, double ds0, int ratio)  throw():
    jmx_(jmx), imx_(imx), jx_(jx), ix_(ix), phic_(phic), lambdac_(lambdac),
    phi1_(phi1), phi2_(phi2), ds0_(ds0), ratio_(ratio),
    pi_(3.14159265358979323846264),
    Earth_radius_(6370997.)
  {
    
  }


  //! Convertion operator.
  /*!
    \param lon longitude.
    \param lat latitude.
    \param j index of the MM5 grid along the East-West direction.
    \param i index of the MM5 grid along the North-South direction.
    \warning Indices order (in MM5) is confusing: j and i are swapped with respect to
    the natural order.
  */
  template <class T>
  void LonlatToMM5LccInd<T>::operator() (const T lon, const T lat,
					 T& j, T& i)
  {

    double ic0, jc0;
    double ic, jc;
    double conv;
    double aux1, aux2, kappa;
    double auxsig;
    double psi1;
    double auxc, yc;
    double x, y, R;
    double auxl;
    double lambdaprima;
    double aux, Rs;
    double auxij;
    double auxtan;

    ic0 = (imx_ + 1.0) / 2.0;
    jc0 = (jmx_ + 1.0) / 2.0;

    ic = (ic0 - ix_) * ratio_ + 0.5;
    jc = (jc0 - jx_) * ratio_ + 0.5;

    conv = 360.0 / (2.0 * pi_);

    aux1 = (45.0 - abs(phi1_) / 2.0) / conv;
    aux2 = (45.0 - abs(phi2_) / 2.0) / conv;
    kappa = ( log10(cos(phi1_/conv)) - log10(cos(phi2_/conv)) )
      / ( log10(tan(aux1)) - log10(tan(aux2)) );

    auxsig = phic_>0 ? 1.0 : -1.0;
    
    psi1 = auxsig * ( pi_ / 2.0 - abs(phi1_) / conv );

    auxc = (auxsig * 90.0 - phic_) / conv / 2.0;
    yc = - (Earth_radius_ / kappa) * sin(psi1)
      * pow(tan(auxc) / tan(psi1/2.0), kappa);

    aux = (auxsig * 90.0 - lat) / (2.0 * conv);
    Rs = (Earth_radius_ / kappa) * sin(psi1) * pow(tan(aux) / tan(psi1/2.0), kappa);

    auxij = kappa * (lon - lambdac_) / conv;
    i = ( ic0 - ( yc/ds0_ + Rs*cos(auxij)/ds0_ ) - ix_ ) * ratio_ - 0.5;

    j = ( jc0 + auxsig * Rs * sin(auxij)/ds0_ - jx_ ) * ratio_ - 0.5;

  }


  //! Default constructor.
  template <class T>
  MM5MercIndToLonlat<T>::MM5MercIndToLonlat(int jmx, int imx, double jx, double ix, double phic,
					    double lambdac, double phi1, double ds, int ratio)  throw():
    jmx_(jmx), imx_(imx), jx_(jx), ix_(ix), phic_(phic), lambdac_(lambdac),
    phi1_(phi1), ds_(ds), ratio_(ratio),
    pi_(3.14159265358979323846264),
    Earth_radius_(6370997.)
  {
    
  }


  //! Convertion operator.
  /*!
    \param j index of the MM5 grid along the East-West direction.
    \param i index of the MM5 grid along the North-South direction.
    \param lon longitude (output).
    \param lat latitude (output).
    \warning Indices order (in MM5) is confusing: j has to be provided first.
  */
  template <class T>
  void MM5MercIndToLonlat<T>::operator() (const T j, const T i,
					  T& lon, T& lat)
  {

    double ic0, jc0;
    double ic, jc;
    double c2;
    double yc;
    double conv;
    double aux;
    double x, y;

    ic0 = (imx_ + 1.0) / 2.0;
    jc0 = (jmx_ + 1.0) / 2.0;

    ic = (ic0 - ix_) * ratio_ + 0.5;
    jc = (jc0 - jx_) * ratio_ + 0.5;

    conv = 180. / pi_;

    c2 = Earth_radius_ * cos(phi1_/conv);
    aux = phic_ / conv;
    yc = c2 * log((1.0 + sin(aux)) / cos(aux));

    aux = (yc + (i + 1 - ic) * ds_) / c2;
    lat = 2.0 * conv * atan(exp(aux)) - 90.0;

    lon = lambdac_ + conv * (j + 1 - jc) * ds_ / c2;

  }


  //! Default constructor.
  template <class T>
  LonlatToMM5MercInd<T>::LonlatToMM5MercInd(int jmx, int imx, double jx, double ix, double phic,
					    double lambdac, double phi1, double ds0, int ratio)  throw():
    jmx_(jmx), imx_(imx), jx_(jx), ix_(ix), phic_(phic), lambdac_(lambdac),
    phi1_(phi1), ds0_(ds0), ratio_(ratio),
    pi_(3.14159265358979323846264),
    Earth_radius_(6370997.)
  {
    
  }


  //! Convertion operator.
  /*!
    \param lon longitude.
    \param lat latitude.
    \param j index of the MM5 grid along the East-West direction.
    \param i index of the MM5 grid along the North-South direction.
    \warning Indices order (in MM5) is confusing: j and i are swapped with respect to
    the natural order.
  */
  template <class T>
  void LonlatToMM5MercInd<T>::operator() (const T lon, const T lat,
					  T& j, T& i)
  {

    double ic0, jc0;
    double ic, jc;
    double c2;
    double yc;
    double conv;
    double aux;
    double x, y;

    ic0 = (imx_ + 1.0) / 2.0;
    jc0 = (jmx_ + 1.0) / 2.0;

    ic = (ic0 - ix_) * ratio_ + 0.5;
    jc = (jc0 - jx_) * ratio_ + 0.5;

    conv = 180. / pi_;

    c2 = Earth_radius_ * cos(phi1_/conv);
    aux = phic_ / conv;
    yc = c2 * log((1.0 + sin(aux)) / cos(aux));

    x = c2 * (lon - lambdac_) / conv;
    aux = lat / conv;
    y = c2 * log((1.0 + sin(aux)) / cos(aux));

    i = (ic0 + (y-yc)/ds0_ - ix_) * ratio_ - 0.5;
    j = (jc0 + x/ds0_ - jx_) * ratio_ - 0.5;

  }


  //! Default constructor.
  template <class T>
  MM5StereIndToLonlat<T>::MM5StereIndToLonlat(int jmx, int imx, double jx, double ix, double phic,
					      double lambdac, double phi1, double ds, int ratio)  throw():
    jmx_(jmx), imx_(imx), jx_(jx), ix_(ix), phic_(phic), lambdac_(lambdac),
    phi1_(phi1), ds_(ds), ratio_(ratio),
    pi_(3.14159265358979323846264),
    Earth_radius_(6370997.)
  {
    
  }


  //! Convertion operator.
  /*!
    \param j index of the MM5 grid along the East-West direction.
    \param i index of the MM5 grid along the North-South direction.
    \param lon longitude (output).
    \param lat latitude (output).
    \warning Indices order (in MM5) is confusing: j has to be provided first.
  */
  template <class T>
  void MM5StereIndToLonlat<T>::operator() (const T j, const T i,
					   T& lon, T& lat)
  {
 
    double ic0, jc0;
    double ic, jc;
    double conv;
    double kappa;
    double auxsig;
    double psi1;
    double auxc, yc;
    double x, y, R;
    double auxl;
    double auxtan;
    double lambdaprima;
    double aux, Rs;
    double auxij;


    ic0 = (imx_ + 1.0) / 2.0;
    jc0 = (jmx_ + 1.0) / 2.0;

    ic = (ic0 - ix_) * ratio_ + 0.5;
    jc = (jc0 - jx_) * ratio_ + 0.5;

    conv = 180. / pi_;

    kappa = 1.0;

    auxsig = phic_<0? -1.0 : 1.0;

    psi1 = auxsig * (pi_/2.0 - abs(phi1_)/conv);

    auxc = (auxsig * 90.0 - phic_) / conv;
    yc = - Earth_radius_ * sin(auxc) * ((1.0 + cos(psi1)) / (1.0 + cos(auxc)));

    x = (jc - j - 1) * ds_;
    y = (i + 1 - ic) * ds_ + yc;
    R = sqrt(x*x + y*y);

    auxl = R/((1.0 + cos(psi1)) * Earth_radius_);
    lat = auxsig * 90.0 - 2.0 * conv * atan(auxl);

    if (y == 0.0)
      {
	auxtan = x<0.? (-pi_/2.) : (pi_/2.);
	lambdaprima = lambdac_ + conv * auxtan / kappa;
      }
      else
	lambdaprima = lambdac_ + conv * atan(x/(auxsig*y)) /kappa;
    
    if (lambdaprima < 180.0)
      lon = lambdaprima + 360.0;
    if ((-180.0 <= lambdaprima) && (lambdaprima <= 180.0))
      lon = lambdaprima;
    if (lambdaprima > 180.0)
      lon = lambdaprima - 360.0;

  }


  //! Default constructor.
  template <class T>
  LonlatToMM5StereInd<T>::LonlatToMM5StereInd(int jmx, int imx, double jx, double ix, double phic,
					      double lambdac, double phi1, double ds0, int ratio)  throw():
    jmx_(jmx), imx_(imx), jx_(jx), ix_(ix), phic_(phic), lambdac_(lambdac),
    phi1_(phi1), ds0_(ds0), ratio_(ratio),
    pi_(3.14159265358979323846264),
    Earth_radius_(6370997.)
  {
    
  }


  //! Convertion operator.
  /*!
    \param lon longitude.
    \param lat latitude.
    \param j index of the MM5 grid along the East-West direction.
    \param i index of the MM5 grid along the North-South direction.
    \warning Indices order (in MM5) is confusing: j and i are swapped with respect to
    the natural order.
  */
  template <class T>
  void LonlatToMM5StereInd<T>::operator() (const T lon, const T lat,
					   T& j, T& i)
  {
 
    double ic0, jc0;
    double ic, jc;
    double conv;
    double kappa;
    double auxsig;
    double psi1;
    double auxc, yc;
    double x, y, R;
    double auxl;
    double auxtan;
    double lambdaprima;
    double aux, Rs;
    double auxij;


    ic0 = (imx_ + 1.0) / 2.0;
    jc0 = (jmx_ + 1.0) / 2.0;

    ic = (ic0 - ix_) * ratio_ + 0.5;
    jc = (jc0 - jx_) * ratio_ + 0.5;

    conv = 180. / pi_;

    kappa = 1.0;

    auxsig = phic_<0? -1.0 : 1.0;

    psi1 = auxsig * (pi_/2.0 - abs(phi1_)/conv);

    auxc = (auxsig * 90.0 - phic_) / conv;
    yc = - Earth_radius_ * sin(auxc) * ((1.0 + cos(psi1)) / (1.0 + cos(auxc)));

    aux = (auxsig * 90.0 - lat) / conv;
    Rs = Earth_radius_ * sin(aux) * ((1.0 + cos(psi1)) / (1.0 + cos(aux)));

    auxij = kappa * (lon - lambdac_) / conv;
    i = (ic0 - (yc / ds0_ + Rs * cos(auxij) / ds0_) - ix_) * ratio_ - 0.5;
    j = (jc0 + auxsig * Rs * sin(auxij) / ds0_ - jx_) * ratio_ - 0.5;

  }


}  // namespace AtmoData.

#define ATMODATA_FILE_COORDTRANSFORM_CXX
#endif
