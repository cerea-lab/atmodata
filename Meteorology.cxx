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


#ifndef ATMODATA_FILE_METEOROLOGY_CXX

#include "Meteorology.hxx"

namespace AtmoData
{


  //! Computes the Richardson number.
  /*!
    Winds may be provided in two ways. The first option is to provide winds on
    interfaces (along x for the zonal wind, along y for the meridional wind).
    The second option is simply to provide winds at nodes (i.e. where the
    potential temperature and the Richardson number are defined).
    \param ZonalWind zonal wind.
    \param MeridionalWind meridional wind.
    \param PotentialTemperature potential temperature.
    \param Richardson (output) Richardson number.
    \param wind_threshold (optional) minimum of the wind shear. Default: 0.1.
  */
  template<class TU, class TV, class TTp, class T, class TG>
  void ComputeRichardson(Data<TU, 4, TG>& ZonalWind,
			 Data<TV, 4, TG>& MeridionalWind,
			 Data<TTp, 4, TG>& PotentialTemperature,
			 Data<T, 4, TG>& Richardson, T wind_threshold)
  {

    int h, i, j, k;

    int Nx = Richardson.GetLength(3);
    int Ny = Richardson.GetLength(2);
    int Nz = Richardson.GetLength(1);
    int Nt = Richardson.GetLength(0);

    Grid<TG>& Levels = Richardson[1];

    const T g(9.81);
    T dudz, dvdz, dwinddz;
    int level;

    if (ZonalWind.GetLength(3) == Nx + 1
	&& MeridionalWind.GetLength(2) == Ny + 1)
      for (h=0; h<Nt; h++)
	for (k=0; k<Nz; k++)
	  for (j=0; j<Ny; j++)
	    for (i=0; i<Nx; i++)
	      {
		
		level = k==0 ? 1 : k;
		
		dudz = 0.5
		  * ( (ZonalWind(h, level, j, i+1)
		       - ZonalWind(h, level-1, j, i+1))
		      / (ZonalWindLevels(h, level, j, i+1)
			 - ZonalWindLevels(h, level-1, j, i+1))
		      + (ZonalWind(h, level, j, i)
			 - ZonalWind(h, level-1, j, i))
		      / (ZonalWindLevels(h, level, j, i)
			 - ZonalWindLevels(h, level-1, j, i)) );
		dvdz = 0.5
		  * ( (MeridionalWind(h, level, j+1, i)
		       - MeridionalWind(h, level-1, j+1, i))
		      / (MeridionalWindLevels(h, level, j+1, i)
			 - MeridionalWindLevels(h, level-1, j+1, i))
		      + (MeridionalWind(h, level, j, i)
			 - MeridionalWind(h, level-1, j, i))
		      / (MeridionalWindLevels(h, level, j, i)
			 - MeridionalWindLevels(h, level-1, j, i)) );
		dwinddz = max(sqrt(dudz*dudz + dvdz*dvdz), wind_threshold);
		Richardson(h, k, j, i) =
		  g * (PotentialTemperature(h, level, j, i)
		       - PotentialTemperature(h, level - 1, j, i))
		  / (dwinddz * dwinddz
		     * PotentialTemperature(h, level - 1, j, i)
		     * (Levels.Value(h, level, j, i)
			- Levels.Value(h, level - 1, j, i)));

	      }
    else
      for (h=0; h<Nt; h++)
	for (k=0; k<Nz; k++)
	  for (j=0; j<Ny; j++)
	    for (i=0; i<Nx; i++)
	      {
		
		level = k==0 ? 1 : k;
		
		dudz = (ZonalWind(h, level, j, i)
			- ZonalWind(h, level-1, j, i))
		  / (ZonalWindLevels(h, level, j, i)
		     - ZonalWindLevels(h, level-1, j, i));
		dvdz = (MeridionalWind(h, level, j, i)
			- MeridionalWind(h, level-1, j, i))
		      / (MeridionalWindLevels(h, level, j, i)
			 - MeridionalWindLevels(h, level-1, j, i));
		dwinddz = max(sqrt(dudz*dudz + dvdz*dvdz), wind_threshold);
		Richardson(h, k, j, i) =
		  g * (PotentialTemperature(h, level, j, i)
		       - PotentialTemperature(h, level - 1, j, i))
		  / (dwinddz * dwinddz
		     * PotentialTemperature(h, level - 1, j, i)
		     * (Levels.Value(h, level, j, i)
			- Levels.Value(h, level - 1, j, i)));

	      }

  }


  //! Computes the surface Richardson number.
  /*!
    Winds may be provided in two ways. The first option is to provide winds on
    interfaces (along x for the zonal wind, along y for the meridional wind).
    The second option is simply to provide winds at nodes (i.e. where the
    potential temperature and the Richardson number are defined).
    \param ZonalWind zonal wind.
    \param MeridionalWind meridional wind.
    \param PotentialTemperature potential temperature.
    \param SurfaceRichardson (output) surface Richardson number.
    \param wind_threshold (optional) minimum of the wind shear. Default: 0.1.
  */
  template<class TU, class TV, class TTp, class T, class TG>
  void ComputeRichardson(Data<TU, 4, TG>& ZonalWind,
			 Data<TV, 4, TG>& MeridionalWind,
			 Data<TTp, 4, TG>& PotentialTemperature,
			 Data<T, 3, TG>& SurfaceRichardson, T wind_threshold)
  {

    int h, i, j;

    int Nx = SurfaceRichardson.GetLength(2);
    int Ny = SurfaceRichardson.GetLength(1);
    int Nt = SurfaceRichardson.GetLength(0);

    Grid<TG>& Levels = PotentialTemperature[1];
    Grid<TG>& MeridionalWindLevels = MeridionalWind[1];
    Grid<TG>& ZonalWindLevels = ZonalWind[1];

    const T g(9.81);
    T dudz, dvdz, dwinddz;

    if (ZonalWind.GetLength(3) == Nx + 1
	&& MeridionalWind.GetLength(2) == Ny + 1)
      for (h=0; h<Nt; h++)
	for (j=0; j<Ny; j++)
	  for (i=0; i<Nx; i++)
	    {
		
	      dudz = 0.5
		* (ZonalWind(h, 0, j, i+1)
		   / ZonalWindLevels.Value(h, 0, j, i+1)
		   + ZonalWind(h, 0, j, i)
		   / ZonalWindLevels.Value(h, 0, j, i));
	      dvdz = 0.5 *
		(MeridionalWind(h, 0, j+1, i)
		 / MeridionalWindLevels.Value(h, 0, j+1, i)
		 + MeridionalWind(h, 0, j, i)
		 / MeridionalWindLevels.Value(h, 0, j, i));
	      dwinddz = max(sqrt(dudz*dudz + dvdz*dvdz), wind_threshold);
	      SurfaceRichardson(h, j, i) =
		g * (PotentialTemperature(h, 1, j, i)
		     - PotentialTemperature(h, 0, j, i))
		/ (dwinddz * dwinddz * PotentialTemperature(h, 0, j, i)
		   * (Levels.Value(h, 1, j, i) - Levels.Value(h, 0, j, i)));
	      
	    }
    else
      for (h=0; h<Nt; h++)
	for (j=0; j<Ny; j++)
	  for (i=0; i<Nx; i++)
	    {
		
	      dudz = ZonalWind(h, 0, j, i)
		/ ZonalWindLevels.Value(h, 0, j, i);
	      dvdz = MeridionalWind(h, 0, j, i)
		/ MeridionalWindLevels.Value(h, 0, j, i);
	      dwinddz = max(sqrt(dudz*dudz + dvdz*dvdz), wind_threshold);
	      SurfaceRichardson(h, j, i) =
		g * (PotentialTemperature(h, 1, j, i)
		     - PotentialTemperature(h, 0, j, i))
		/ (dwinddz * dwinddz * PotentialTemperature(h, 0, j, i)
		   * (Levels.Value(h, 1, j, i) - Levels.Value(h, 0, j, i)));

	    }

  }


  //! Computes the potential temperature.
  /*!
    Formula: PotentialTemperature = Temperature * (Pressure / P0)^(-r/cp).
    \param Temperature temperature (or virtual temperature).
    \param Pressure pressure.
    \param PotentialTemperature (output) potential temperature.
    \param P0 (optional) standard pressure. Default: 101325 Pa.
    \param cp (optional) specific heat of dry air at constant pressure.
    Default: 1005 J.kg^{-1}.K^{-1}.
    \param r (optional) molar gas constant for air.
    Default: 287.0 J.kg^{-1}.K^{-1}.
  */
  template<class TT, class TP, class T, class TG>
  void ComputePotentialTemperature(Data<TT, 4, TG>& Temperature,
				   Data<TP, 4, TG>& Pressure,
				   Data<T, 4, TG>& PotentialTemperature,
				   T P0, T cp, T r)
  {

    int h, i, j, k;

    int Nx = PotentialTemperature.GetLength(3);
    int Ny = PotentialTemperature.GetLength(2);
    int Nz = PotentialTemperature.GetLength(1);
    int Nt = PotentialTemperature.GetLength(0);

    T ratio = -r / cp;

    for (h=0; h<Nt; h++)
      for (k=0; k<Nz; k++)
	for (j=0; j<Ny; j++)
	  for (i=0; i<Nx; i++)
	    PotentialTemperature(h, k, j, i) =
	      Temperature(h, k, j, i)
	      * pow(Pressure(h, k, j, i) / P0, ratio);
    
  }


  //! Computes the potential temperature.
  /*!
    Formula: PotentialTemperature = Temperature * (Pressure / P0)^(-r/cp).
    \param Temperature temperature (or virtual temperature).
    \param Pressure pressure.
    \param PotentialTemperature (output) potential temperature.
    \param P0 (optional) standard pressure. Default: 101325 Pa.
    \param cp (optional) specific heat of dry air at constant pressure.
    Default: 1005 J.kg^{-1}.K^{-1}.
    \param r (optional) molar gas constant for air.
    Default: 287.0 J.kg^{-1}.K^{-1}.
  */
  template<class TT, class TP, class T, class TG>
  void ComputePotentialTemperature(Data<TT, 3, TG>& Temperature,
				   Data<TP, 3, TG>& Pressure,
				   Data<T, 3, TG>& PotentialTemperature,
				   T P0, T cp, T r)
  {

    int h, i, j;

    int Nx = PotentialTemperature.GetLength(2);
    int Ny = PotentialTemperature.GetLength(1);
    int Nt = PotentialTemperature.GetLength(0);

    T ratio = -r / cp;

    for (h=0; h<Nt; h++)
      for (j=0; j<Ny; j++)
	for (i=0; i<Nx; i++)
	  PotentialTemperature(h, j, i) =
	    Temperature(h, j, i) * pow(Pressure(h, j, i) / P0, ratio);
    
  }


  //! Computes the relative humidity from the specific humidity.
  /*!
    \param SpecificHumidity specific humidity (kg/kg).
    \param Temperature temperature (or virtual temperature) (K).
    \param Pressure pressure (Pa).
    \param RelativeHumidity (output) relative humidity.
  */
  template<class TS, class TT, class TP, class T, class TG>
  void ComputeRelativeHumidity(Data<TS, 4, TG>& SpecificHumidity,
			       Data<TT, 4, TG>& Temperature,
			       Data<TP, 4, TG>& Pressure,
			       Data<T, 4, TG>& RelativeHumidity)
  {
    int h, k, j, i;
    int Nt(RelativeHumidity.GetLength(0));
    int Nz(RelativeHumidity.GetLength(1));
    int Ny(RelativeHumidity.GetLength(2));
    int Nx(RelativeHumidity.GetLength(3));

    T P_sat;
    for (h = 0; h < Nt; h++)
      for (k = 0; k < Nz; k++)
	for (j = 0; j < Ny; j++)
	  for (i = 0; i < Nx; i++)
	    {
	      P_sat = 611.2 * exp(17.67 * (Temperature(h, k, j, i) - 273.15)
				  / (Temperature(h, k, j, i) - 29.65));
	      RelativeHumidity(h, k, j, i) = SpecificHumidity(h, k, j, i)
		* Pressure(h, k, j, i)
		/ ( (0.62197 * (1.0 - SpecificHumidity(h, k, j, i))
		     + SpecificHumidity(h, k, j, i) ) * P_sat);
	    }
  }


  //! Computes the critival relative humidity.
  /*!
    Formula: CriticalRelativeHumidity = 1.0 - coeff0 * sig
    * (1.0 - sig) * (1.0 + (sig - 0.5) * coeff1) where
    sig = Pressure / SurfacePressure.
    \param SurfacePressure surface pressure (Pa).
    \param Pressure pressure (Pa).
    \param CriticalRelativeHumidity (output) critical relative humidity.
    \param coeff0 coefficient (see the formula). Default: 2.0.
    \param coeff1 coefficient (see the formula). Default: sqrt(3.).
  */
  template<class TS, class TP, class T, class TG>
  void ComputeCriticalRelativeHumidity(Data<TS, 3, TG>& SurfacePressure,
				       Data<TP, 4, TG>& Pressure,
				       Data<T, 4, TG>& CriticalRelativeHumidity,
				       T coeff0, T coeff1)
  {
    int h, k, j, i;
    int Nt(CriticalRelativeHumidity.GetLength(0));
    int Nz(CriticalRelativeHumidity.GetLength(1));
    int Ny(CriticalRelativeHumidity.GetLength(2));
    int Nx(CriticalRelativeHumidity.GetLength(3));

    T sig;
    for (h = 0; h < Nt; h++)
      for (k = 0; k < Nz; k++)
	for (j = 0; j < Ny; j++)
	  for (i = 0; i < Nx; i++)
	    {
	      sig = Pressure(h, k, j, i) / SurfacePressure(h, j, i);
	      CriticalRelativeHumidity(h, k, j, i) = 1.0 - coeff0 * sig
		* (1.0 - sig) * (1.0 + (sig - 0.5) * coeff1);
	    }
  }


  //! Computes the critival relative humidity.
  /*!
    Extended formula: CriticalRelativeHumidity = 1.0 - coeff0 * sig^a0
    * (1.0 - sig)^a1 * (1.0 + (sig - 0.5) * coeff1) where
    sig = Pressure / SurfacePressure.
    \param SurfacePressure surface pressure (Pa).
    \param Pressure pressure (Pa).
    \param CriticalRelativeHumidity (output) critical relative humidity.
    \param coeff0 coefficient (see the formula). Default: 1.1.
    \param coeff1 coefficient (see the formula). Default: sqrt(1.3).
    \param a0 exponent (see the formula). Default: 0.0.
    \param a1 exponent (see the formula). Default: 1.1.
  */
  template<class TS, class TP, class T, class TG>
  void ComputeCriticalRelativeHumidity_extended(Data<TS, 3, TG>& SurfacePressure,
						Data<TP, 4, TG>& Pressure,
						Data<T, 4, TG>& CriticalRelativeHumidity,
						T coeff0, T coeff1,
						T a0, T a1)
  {
    int h, k, j, i;
    int Nt(CriticalRelativeHumidity.GetLength(0));
    int Nz(CriticalRelativeHumidity.GetLength(1));
    int Ny(CriticalRelativeHumidity.GetLength(2));
    int Nx(CriticalRelativeHumidity.GetLength(3));

    T sig;
    for (h = 0; h < Nt; h++)
      for (k = 0; k < Nz; k++)
	for (j = 0; j < Ny; j++)
	  for (i = 0; i < Nx; i++)
	    {
	      sig = Pressure(h, k, j, i) / SurfacePressure(h, j, i);
	      CriticalRelativeHumidity(h, k, j, i) =
		1.0 - coeff0 * pow(sig, a0)
		* pow(T(1.) - sig, a1) * (1.0 + (sig - 0.5) * coeff1);
	    }
  }


  //! Computes the critical relative humidity.
  /*!
    Formula: inside the boundary layer,
    CriticalRelativeHumidity = BL_CRH and, above the boundary layer,
    CriticalRelativeHumidity = 1.0 - coeff0 * sig
    * (1.0 - sig) * (1.0 + (sig - 0.5) * coeff1) where
    sig = Pressure / SurfacePressure.
    \param BoundaryLayerHeight boundary layer height (m).
    \param SurfacePressure surface pressure (Pa).
    \param Pressure pressure (Pa).
    \param CriticalRelativeHumidity (output) critical relative humidity.
    \param coeff0 coefficient (see the formula). Default: 2.0.
    \param coeff1 coefficient (see the formula). Default: sqrt(3.).
    \param BL_CRH critical relative humidity within the boundary layer.
    Default: 0.98.
  */
  template<class TB, class TS, class TP, class T, class TG>
  void ComputeCriticalRelativeHumidity(Data<TB, 3, TG>& BoundaryLayerHeight,
				       Data<TS, 3, TG>& SurfacePressure,
				       Data<TP, 4, TG>& Pressure,
				       Data<T, 4, TG>& CriticalRelativeHumidity,
				       T coeff0, T coeff1, T BL_CRH)
  {
    int h, k, j, i;
    int Nt(CriticalRelativeHumidity.GetLength(0));
    int Nz(CriticalRelativeHumidity.GetLength(1));
    int Ny(CriticalRelativeHumidity.GetLength(2));
    int Nx(CriticalRelativeHumidity.GetLength(3));

    T sig;
    for (h = 0; h < Nt; h++)
      for (k = 0; k < Nz; k++)
	for (j = 0; j < Ny; j++)
	  for (i = 0; i < Nx; i++)
	    if (CriticalRelativeHumidity[1].Value(h, k, j, i)
		< BoundaryLayerHeight(h, j, i))
	      CriticalRelativeHumidity(h, k, j, i) = BL_CRH;
	    else
	      {
		sig = Pressure(h, k, j, i) / SurfacePressure(h, j, i);
		CriticalRelativeHumidity(h, k, j, i) = 1.0 - coeff0 * sig
		  * (1.0 - sig) * (1.0 + (sig - 0.5) * coeff1);
	      }
  }


  //! Computes the critical relative humidity.
  /*!
    The relative humidity is set to CRH_0 if Pressure > P_0,
    to CRH_1 if P_0 >= Pressure > P_1 and to CRH_2 otherwise.
    \param Pressure pressure (Pa).
    \param CriticalRelativeHumidity (output) critical relative humidity.
    \param CRH_0 critical relative humidity in the first layer. Default: 0.75.
    \param CRH_1 critical relative humidity in the second layer.
    Default: 0.95.
    \param CRH_2 critical relative humidity in the third layer. Default: 0.95.
    \param P_0 first pressure limit. Default: 70 000 Pa.
    \param P_1 second pressure limit. Default: 40 000 Pa.
  */
  template<class TP, class T, class TG>
  void ComputeCriticalRelativeHumidity(Data<TP, 4, TG>& Pressure,
				       Data<T, 4, TG>& CriticalRelativeHumidity,
				       T CRH_0, T CRH_1, T CRH_2,
				       T P_0, T P_1)
  {
    int h, k, j, i;
    int Nt(CriticalRelativeHumidity.GetLength(0));
    int Nz(CriticalRelativeHumidity.GetLength(1));
    int Ny(CriticalRelativeHumidity.GetLength(2));
    int Nx(CriticalRelativeHumidity.GetLength(3));

    for (h = 0; h < Nt; h++)
      for (k = 0; k < Nz; k++)
	for (j = 0; j < Ny; j++)
	  for (i = 0; i < Nx; i++)
	    if (Pressure(h, k, j, i) > P_0)
	      CriticalRelativeHumidity(h, k, j, i) = CRH_0;
	    else if (Pressure(h, k, j, i) > P_1)
	      CriticalRelativeHumidity(h, k, j, i) = CRH_1;
	    else
	      CriticalRelativeHumidity(h, k, j, i) = CRH_2;
  }


  //! Computes the cloud fraction.
  /*!
    Formula: CloudFraction = [ (RelativeHumidity - CriticalRelativeHumidity)
    / (1.0 - CriticalRelativeHumidity) ]^2
    \param RelativeHumidity relative humidity.
    \param CriticalRelativeHumidity critical relative humidity.
    \param CloudFraction (output) cloud fraction.
  */
  template<class TR, class TC, class T, class TG>
  void ComputeCloudFraction(Data<TR, 4, TG>& RelativeHumidity,
			    Data<TC, 4, TG>& CriticalRelativeHumidity,
			    Data<T, 4, TG>& CloudFraction)
  {
    int h, k, j, i;
    int Nt(CloudFraction.GetLength(0));
    int Nz(CloudFraction.GetLength(1));
    int Ny(CloudFraction.GetLength(2));
    int Nx(CloudFraction.GetLength(3));

    T tmp, crh;
    for (h = 0; h < Nt; h++)
      for (k = 0; k < Nz; k++)
	for (j = 0; j < Ny; j++)
	  for (i = 0; i < Nx; i++)
	    {
	      crh = CriticalRelativeHumidity(h, k, j, i);
	      tmp = RelativeHumidity(h, k, j, i) - crh;
	      if (tmp < 0)
		CloudFraction(h, k, j, i) = 0.;
	      else
		{
		  tmp = tmp / (1. - crh);
		  CloudFraction(h, k, j, i) = tmp * tmp;
		}
	    }
  }


  //! Computes the cloud fraction.
  /*!
    Formula: inside the boundary layer,
    CloudFraction = 0.34 * [ (RelativeHumidity - CriticalRelativeHumidity)
    / (1.0 - CriticalRelativeHumidity) ] and, above the boundary layer,
    CloudFraction = [ (RelativeHumidity - CriticalRelativeHumidity)
    / (1.0 - CriticalRelativeHumidity) ]^2
    \param BoundaryLayerHeight boundary layer height (m).
    \param RelativeHumidity relative humidity.
    \param CriticalRelativeHumidity critical relative humidity.
    \param CloudFraction (output) cloud fraction.
  */
  template<class TP, class TR, class TC, class T, class TG>
  void ComputeCloudFraction(Data<TP, 3, TG>& BoundaryLayerHeight,
			    Data<TR, 4, TG>& RelativeHumidity,
			    Data<TC, 4, TG>& CriticalRelativeHumidity,
			    Data<T, 4, TG>& CloudFraction)
  {
    int h, k, j, i;
    int Nt(CloudFraction.GetLength(0));
    int Nz(CloudFraction.GetLength(1));
    int Ny(CloudFraction.GetLength(2));
    int Nx(CloudFraction.GetLength(3));

    T tmp, crh;
    for (h = 0; h < Nt; h++)
      for (k = 0; k < Nz; k++)
	for (j = 0; j < Ny; j++)
	  for (i = 0; i < Nx; i++)
	    {
	      crh = CriticalRelativeHumidity(h, k, j, i);
	      tmp = RelativeHumidity(h, k, j, i) - crh;
	      if (tmp < 0)
		CloudFraction(h, k, j, i) = 0.;
	      else
		if (CloudFraction[1].Value(h, k, j, i)
		    < BoundaryLayerHeight(h, j, i))
		  CloudFraction(h, k, j, i) = 0.34 * tmp / (1. - crh);
		else
		  {
		    tmp = tmp / (1. - crh);
		    CloudFraction(h, k, j, i) = tmp * tmp;
		  }
	    }
  }


  //! Computes the module of a 2D-vectors field.
  /*!
    This function was initially dedicated to winds. In this case, zonal winds
    and meridional winds are provided and the module of the wind is computed
    (assuming that the vertical wind is zero). Winds may be provided in two
    ways. The first option is to provide winds on interfaces (along x for
    the zonal wind, along y for the meridional wind). The second option is
    simply to provide winds at nodes (i.e. where the module is defined).
    \param U first component of vectors.
    \param V second component of vectors.
    \param Module (output) module.
  */
  template<class TU, class TV, class T, class TG>
  void ComputeModule(Data<TU, 4, TG>& U, Data<TV, 4, TG>& V,
		     Data<T, 4, TG>& Module)
  {

    int h, i, j, k;

    int Nx = Module.GetLength(3);
    int Ny = Module.GetLength(2);
    int Nz = Module.GetLength(1);
    int Nt = Module.GetLength(0);

    T u, v;

    if ( (U.GetLength(3) == Nx + 1) && (V.GetLength(2) == Ny + 1) )
      for (h=0; h<Nt; h++)
	for (k=0; k<Nz; k++)
	  for (j=0; j<Ny; j++)
	    for (i=0; i<Nx; i++)
	      {
		u = 0.5 * (U(h, k, j, i+1) + U(h, k, j, i));
		v = 0.5 * (V(h, k, j+1, i) + V(h, k, j, i));
		Module(h, k, j, i) = sqrt(u*u + v*v);
	      }
    else
      for (h=0; h<Nt; h++)
	for (k=0; k<Nz; k++)
	  for (j=0; j<Ny; j++)
	    for (i=0; i<Nx; i++)
	      {
		u = U(h, k, j, i);
		v = V(h, k, j, i);
		Module(h, k, j, i) = sqrt(u*u + v*v);
	      }

  }


  //! Computes the module of a 2D-vectors field on the surface.
  /*!
    This function was initially dedicated to winds. In this case, zonal winds
    and meridional winds are provided and the module of the wind is computed
    (assuming that the vertical wind is zero). Winds may be provided in two
    ways. The first option is to provide winds on interfaces (along x
    for the zonal wind, along y for the meridional wind). The second option is
    simply to provide winds at nodes (i.e. where the module is defined).
    \param U first component of vectors.
    \param V second component of vectors.
    \param Module (output) surface module.
  */
  template<class TU, class TV, class T, class TG>
  void ComputeModule(Data<TU, 4, TG>& U, Data<TV, 4, TG>& V,
		     Data<T, 3, TG>& Module)
  {

    int h, i, j;

    int Nx = Module.GetLength(2);
    int Ny = Module.GetLength(1);
    int Nt = Module.GetLength(0);

    T u, v;

    if ( (U.GetLength(3) == Nx + 1) && (V.GetLength(2) == Ny + 1) )
      for (h=0; h<Nt; h++)
	for (j=0; j<Ny; j++)
	  for (i=0; i<Nx; i++)
	    {
	      u = 0.5 * (U(h, 0, j, i+1) + U(h, 0, j, i));
	      v = 0.5 * (V(h, 0, j+1, i) + V(h, 0, j, i));
	      Module(h, j, i) = sqrt(u*u + v*v);
	    }
    else
      for (h=0; h<Nt; h++)
	for (j=0; j<Ny; j++)
	  for (i=0; i<Nx; i++)
	    {
	      u = U(h, 0, j, i);
	      v = V(h, 0, j, i);
	      Module(h, j, i) = sqrt(u*u + v*v);
	    }

  }


  //! Computes the module of a 2D-vectors field on the surface.
  /*!
    This function was initially dedicated to winds. In this case, zonal winds
    and meridional winds are provided and the module of the wind is computed
    (assuming that the vertical wind is zero). Winds may be provided in two
    ways. The first option is to provide winds on interfaces (along x for
    the zonal wind, along y for the meridional wind). The second option is
    simply to provide winds at nodes (i.e. where the module is defined).
    \param U first component of vectors.
    \param V second component of vectors.
    \param Module (output) surface module.
  */
  template<class TU, class TV, class T, class TG>
  void ComputeModule(Data<TU, 3, TG>& U, Data<TV, 3, TG>& V,
		     Data<T, 3, TG>& Module)
  {

    int h, i, j;

    int Nx = Module.GetLength(2);
    int Ny = Module.GetLength(1);
    int Nt = Module.GetLength(0);

    T u, v;

    if ( (U.GetLength(2) == Nx + 1) && (V.GetLength(1) == Ny + 1) )
      for (h=0; h<Nt; h++)
	for (j=0; j<Ny; j++)
	  for (i=0; i<Nx; i++)
	    {
	      u = 0.5 * (U(h, j, i+1) + U(h, j, i));
	      v = 0.5 * (V(h, j+1, i) + V(h, j, i));
	      Module(h, j, i) = sqrt(u*u + v*v);
	    }
    else
      for (h=0; h<Nt; h++)
	for (j=0; j<Ny; j++)
	  for (i=0; i<Nx; i++)
	    {
	      u = U(h, j, i);
	      v = V(h, j, i);
	      Module(h, j, i) = sqrt(u*u + v*v);
	    }

  }


  //! Computes low, medium and high cloudiness, cloud base and top.
  /*!
    \param CloudFraction cloud fraction.
    \param Pressure pressure (Pa).
    \param GridZ_interf altitudes of interfaces (m).
    \param LowIndices vertical indices of base and top of low clouds.
    \param MediumIndices vertical indices of base and top of medium clouds.
    \param HighIndices vertical indices of base and top of high clouds.
    \param LowCloudiness low cloudiness.
    \param MediumCloudiness medium cloudiness.
    \param HighCloudiness high cloudiness.
    \param P_0 first pressure limit. Default: 80 000 Pa.
    \param P_1 second pressure limit. Default: 45 000 Pa.
  */
  template<class TC, class TP, class T, class TG>
  void ComputeCloudiness(Data<TC, 4, TG>& CloudFraction,
			 Data<TP, 4, TG>& Pressure,
			 Grid<TG>& GridZ_interf,
			 Data<int, 4>& LowIndices,
			 Data<int, 4>& MediumIndices,
			 Data<int, 4>& HighIndices,
			 Data<T, 3, TG>& LowCloudiness,
			 Data<T, 3, TG>& MediumCloudiness,
			 Data<T, 3, TG>& HighCloudiness,
			 T P_0, T P_1)
  {
    int h, k, j, i;
    int Nt(CloudFraction.GetLength(0));
    int Nz(CloudFraction.GetLength(1));
    int Ny(CloudFraction.GetLength(2));
    int Nx(CloudFraction.GetLength(3));

    LowCloudiness.SetZero();
    MediumCloudiness.SetZero();
    HighCloudiness.SetZero();

    TC cloud_max;
    bool above, below;
    int k_base, k_top, k_max;
    for (h = 0; h < Nt; h++)
      for (j = 0; j < Ny; j++)
	for (i = 0; i < Nx; i++)
	  {

	    /*** Low clouds ***/

	    cloud_max = 0;
	    // The first level is excluded.
	    for (k = 1; k < Nz && Pressure(h, k, j, i) > P_0; k++)
	      cloud_max = max(cloud_max, CloudFraction(h, k, j, i));
	    below = true; above = false;
	    k_base = 0; k_top = 0;
	    for (k = 1; k < Nz && Pressure(h, k, j, i) > P_0 && !above; k++)
	      {
		below = below && ( CloudFraction(h, k, j, i) < 0.5 * cloud_max
				   || CloudFraction(h, k, j, i) == 0 );
		above = !below && CloudFraction(h, k, j, i) < 0.5 * cloud_max;
		if (!below && k_base == 0)
		  k_base = k;
		if (above)
		  k_top = k;
	      }
	    if (above)
	      while (k < Nz && Pressure(h, k, j, i) > P_0)
		k++;
	    k_max = k - 1;
	    // Goes up to P_0.
	    if (k_base > k_top)
	      k_top = k_max + 1;
	    LowIndices(h, j, i, 0) = k_base;
	    LowIndices(h, j, i, 1) = k_top;
	    k = k_base;
	    // k_top == 0 means no cloud.
	    while (k < k_top && k_top != 0)
	      {
		LowCloudiness(h, j, i) += CloudFraction(h, k, j, i)
		  * (GridZ_interf.Value(h, k + 1, j, i)
		     - GridZ_interf.Value(h, k, j, i));
		k++;
	      }
	    if (k_top != 0)
	      LowCloudiness(h, j, i) /=
		GridZ_interf.Value(h, k_top, j, i)
		- GridZ_interf.Value(h, k_base, j, i);

	    /*** Medium clouds ***/

	    cloud_max = 0;
	    // Starts above low clouds.
	    for (k = k_max + 1; k < Nz && Pressure(h, k, j, i) > P_1; k++)
	      cloud_max = max(cloud_max, CloudFraction(h, k, j, i));
	    below = true; above = false;
	    k_base = 0; k_top = 0;
	    for (k = k_max + 1; k < Nz && Pressure(h, k, j, i) > P_1
		   && !above; k++)
	      {
		below = below && ( CloudFraction(h, k, j, i) < 0.5 * cloud_max
				   || CloudFraction(h, k, j, i) == 0 );
		above = !below && CloudFraction(h, k, j, i) < 0.5 * cloud_max;
		if (!below && k_base == 0)
		  k_base = k;
		if (above)
		  k_top = k;
	      }
	    if (above)
	      while (k < Nz && Pressure(h, k, j, i) > P_1)
		k++;
	    k_max = k - 1;
	    // Goes up to P_1.
	    if (k_base > k_top)
	      k_top = k_max + 1;
	    MediumIndices(h, j, i, 0) = k_base;
	    MediumIndices(h, j, i, 1) = k_top;
	    k = k_base;
	    // k_top == 0 means no cloud.
	    while (k < k_top && k_top != 0)
	      {
		MediumCloudiness(h, j, i) += CloudFraction(h, k, j, i)
		  * (GridZ_interf.Value(h, k + 1, j, i)
		     - GridZ_interf.Value(h, k, j, i));
		k++;
	      }
	    if (k_top != 0)
	      MediumCloudiness(h, j, i) /=
		GridZ_interf.Value(h, k_top, j, i)
		- GridZ_interf.Value(h, k_base, j, i);

	    /*** High clouds ***/

	    cloud_max = 0;
	    // Starts above low clouds.
	    for (k = k_max + 1; k < Nz; k++)
	      cloud_max = max(cloud_max, CloudFraction(h, k, j, i));
	    below = true; above = false;
	    k_base = 0; k_top = 0;
	    for (k = k_max + 1; k < Nz && !above; k++)
	      {
		below = below && ( CloudFraction(h, k, j, i) < 0.5 * cloud_max
				   || CloudFraction(h, k, j, i) == 0 );
		above = !below && CloudFraction(h, k, j, i) < 0.5 * cloud_max;
		if (!below && k_base == 0)
		  k_base = k;
		if (above)
		  k_top = k;
	      }
	    k_max = k - 1;
	    // Goes up to the top.
	    if (k_base > k_top)
	      k_top = k_max + 1;
	    HighIndices(h, j, i, 0) = k_base;
	    HighIndices(h, j, i, 1) = k_top;
	    k = k_base;
	    // k_top == 0 means no cloud.
	    while (k < k_top && k_top != 0)
	      {
		HighCloudiness(h, j, i) += CloudFraction(h, k, j, i)
		  * (GridZ_interf.Value(h, k + 1, j, i)
		     - GridZ_interf.Value(h, k, j, i));
		k++;
	      }
	    if (k_top != 0)
	      HighCloudiness(h, j, i) /=
		GridZ_interf.Value(h, k_top, j, i)
		- GridZ_interf.Value(h, k_base, j, i);
	  }
  }


  //! Computes the height of cloud basis.
  /*!
    \param Temperature temperature (K).
    \param Pressure pressure (Pa).
    \param Humidity relative humidity (kg/kg).
    \param CriticalRelativeHumidity function that returns the critical
    relative humidity as function of the altitude, the pressure
    and reference pressure.
    \param CloudHeight (output) altitudes of cloud basis.
  */
  template <class TP, class TH,
	    class T, class TG>
  void ComputeCloudHeight(Data<TP, 4, TG>& Pressure,
			  Data<TH, 4, TG>& Humidity,
			  T (CriticalRelativeHumidity)(const T&, const T&,
						       const T&),
			  Data<T, 3, TG>& CloudHeight)
  {

    int h, k, j, i;
    int Nt(CloudHeight.GetLength(0));
    int Nz(Pressure.GetLength(1));
    int Ny(CloudHeight.GetLength(1));
    int Nx(CloudHeight.GetLength(2));    

    // Index "0" and "1" refer to two contiguous levels.
    T rh0, rh1, rhc;

    T max_height = 2. * Pressure[1].Value(0, Nz-1, 0, 0);
    CloudHeight.Fill(max_height);

    for (h=0; h<Nt; h++)
      for (j=0; j<Ny; j++)
	for (i=0; i<Nx; i++)
	  {

	    rh0 = Humidity(h, Nz-1, j, i);
	    
	    k = 0;
	    while ( (k<Nz) && (CloudHeight(h, j, i) == max_height) )
	      {

		rh1 = Humidity(h, k, j, i);

		// Critical relative humidity.
		rhc = CriticalRelativeHumidity(Pressure[1].Value(h, k, j, i),
					       Pressure(h, k, j, i),
					       Pressure(h, 0, j, i));

		if (rh1 >= rhc)  // Above a cloud.
		  CloudHeight(h, j, i) = Pressure[1].Value(h, k, j, i);

		// For the next level.
		rh0 = rh1;

		k++;

	      }

	  }
  }


  //! Computes the height of cloud basis.
  /*!
    \param Humidity relative humidity (kg/kg).
    \param CriticalRelativeHumidity critical relative humidity.
    \param CloudHeight (output) altitudes of cloud basis.
  */
  template <class TH, class TCRH, class T, class TG>
  void ComputeCloudHeight(Data<TH, 4, TG>& Humidity,
			  Data<TCRH, 4, TG>& CriticalRelativeHumidity,
			  Data<T, 3, TG>& CloudHeight)
  {

    int h, k, j, i;
    int Nt(CloudHeight.GetLength(0));
    int Nz(Humidity.GetLength(1));
    int Ny(CloudHeight.GetLength(1));
    int Nx(CloudHeight.GetLength(2));    

    // Index "0" and "1" refer to two contiguous levels.
    T rh0, rh1, rhc;

    T max_height = 2. * Humidity[1].Value(0, Nz-1, 0, 0);
    CloudHeight.Fill(max_height);

    for (h=0; h<Nt; h++)
      for (j=0; j<Ny; j++)
	for (i=0; i<Nx; i++)
	  {

	    rh0 = Humidity(h, Nz-1, j, i);
	    
	    k = 0;
	    while ( (k<Nz) && (CloudHeight(h, j, i) == max_height) )
	      {

		rh1 = Humidity(h, k, j, i);

		// Critical relative humidity.
		rhc = CriticalRelativeHumidity(h, k, j, i);

		if (rh1 >= rhc)  // Above a cloud.
		  CloudHeight(h, j, i) = Humidity[1].Value(h, k, j, i);

		// For the next level.
		rh0 = rh1;

		k++;

	      }

	  }
  }


  //! Computes the height of cloud basis.
  /*!
    \param LowIndices vertical indices of base and top of low clouds.
    \param MediumIndices vertical indices of base and top of medium clouds.
    \param HighIndices vertical indices of base and top of high clouds.
    \param CloudHeight (output) altitudes of cloud basis.
  */
  template <class T, class TG>
  void ComputeCloudHeight(Data<int, 4>& LowIndices,
			  Data<int, 4>& MediumIndices,
			  Data<int, 4>& HighIndices,
			  Grid<TG>& GridZ_interf,
			  Data<T, 3, TG>& CloudHeight)
  {
    int h, j, i;
    int Nt(CloudHeight.GetLength(0));
    int Ny(CloudHeight.GetLength(1));
    int Nx(CloudHeight.GetLength(2));

    CloudHeight.SetZero();

    for (h = 0; h < Nt; h++)
      for (j = 0; j < Ny; j++)
	for (i = 0; i < Nx; i++)
	  if (LowIndices(h, j, i, 0) != 0)
	    CloudHeight(h, j, i) =
	      GridZ_interf.Value(h, LowIndices(h, j, i, 0), j, i);
	  else if (MediumIndices(h, j, i, 0) != 0)
	    CloudHeight(h, j, i) =
	      GridZ_interf.Value(h, MediumIndices(h, j, i, 0), j, i);
	  else if (HighIndices(h, j, i, 0) != 0)
	    CloudHeight(h, j, i) =
	      GridZ_interf.Value(h, HighIndices(h, j, i, 0), j, i);
  }


  //! Computes the pressure from the surface pressure.
  /*!
    Formula: Pressure_k = alpha_k * P0 + beta_k * SurfacePressure,
    where k is the level index.
    \param alpha coefficients.
    \param beta coefficients.
    \param SurfacePressure surface pressure.
    \param Pressure (output) pressure.
    \param P0 (optional) standard pressure. Default: 101325 Pa.
  */
  template <class Ta, class Tb, class TSP,
	    class T, class TG>
  void ComputePressure(Data<Ta, 1, TG>& alpha, Data<Tb, 1, TG>& beta,
		       Data<TSP, 3, TG>& SurfacePressure,
		       Data<T, 4, TG>& Pressure, T P0)
  {

    int h, i, j, k;

    int Nx = Pressure.GetLength(3);
    int Ny = Pressure.GetLength(2);
    int Nz = Pressure.GetLength(1);
    int Nt = Pressure.GetLength(0);

    for (h=0; h<Nt; h++)
      for (k=0; k<Nz; k++)
	for (j=0; j<Ny; j++)
	  for (i=0; i<Nx; i++)
	    Pressure(h, k, j, i) = alpha(k) * P0
	      + beta(k) * SurfacePressure(h, j, i);

  }


  //! Computes the altitudes from pressure and temperature fields.
  /*!
    Level heights are computed according to:
    Z_{k+1} = Z_k + (r * T_k / g) * log(P_k/P_{k+1})
    where Z is the altitude, T the temperature, P the pressure,
    k the level index, r the molar gas constant for dry air
    and g the standard gravity.
    \par For the first level, Z_0 = r * T_0 / g * log(PS / P_0)
    where PS is the surface pressure.
    \param SurfacePressure surface pressure (Pa).
    \param Pressure pressure (Pa).
    \param Height (output) altitudes (m).
    \param g (optional) standard gravity. Default: 9.80665.
    \param r (optional) molar gas constant for dry air. Default: 287.0.
    \note Temperature, Pressure and Height must be defined on the same grid.
  */
  template<class TPS, class TP, class TT, class T, class TG>
  void ComputeHeight(Data<TPS, 3, TG>& SurfacePressure,
		     Data<TP, 4, TG>& Pressure,
		     Data<TT, 4, TG>& Temperature,
		     Grid<T>& Height, T g, T r)
  {

    int h, i, j, k;

    int Nx = Height.GetLength(3);
    int Ny = Height.GetLength(2);
    int Nz = Height.GetLength(1);
    int Nt = Height.GetLength(0);

    for (h=0; h<Nt; h++)
      for (k=0; k<Nz; k++)
	for (j=0; j<Ny; j++)
	  for (i=0; i<Nx; i++)
	    if (k==0)
	      Height.Value(h, k, j, i) = r / g * Temperature(h, k, j, i)
		* log(SurfacePressure(h, j, i) / Pressure(h, k, j, i));
	    else
	      Height.Value(h, k, j, i) = Height.Value(h, k-1, j, i)
		- r / g * Temperature(h, k, j, i)
		* log(Pressure(h, k, j, i) / Pressure(h, k-1, j, i));

  }


  //! Computes the altitudes at interfaces from pressure and
  //! temperature fields.
  /*!
    Level heights are computed according to:
    Z_{k+1/2} = Z_{k-1/2} - (r * T_k / g) * log(P_{k+1/2}/P_{k-1/2})
    where Z is the altitude, T the temperature, P the pressure,
    k the level index, r the molar gas constant for dry air
    and g the standard gravity.
    \param Pressure pressure (Pa).
    \param Height (output) altitudes (m).
    \param g (optional) standard gravity. Default: 9.80665.
    \param r (optional) molar gas constant for dry air. Default: 287.0.
    \param ground_set (optional) true if ground-level altitudes are set,
    false if they have to be set (they are set to zero in this case).
    Default: false.
    \note Temperature is provided at middle points (not interfaces).
    Pressure and Height are defined at interfaces (including ground-level).
  */
  template<class TP, class TT, class T, class TG>
  void ComputeInterfHeight(Data<TP, 4, TG>& Pressure,
			   Data<TT, 4, TG>& Temperature,
			   Grid<T>& Height, bool ground_set, T g, T r)
  {

    int h, i, j, k;

    int Nx = Height.GetLength(3);
    int Ny = Height.GetLength(2);
    int Nz = Height.GetLength(1) - 1;
    int Nt = Height.GetLength(0);

    if (!ground_set)
      for (h=0; h<Nt; h++)
	for (j=0; j<Ny; j++)
	  for (i=0; i<Nx; i++)
	    Height.Value(h, 0, j, i) = T(0);

    for (h=0; h<Nt; h++)
      for (k=0; k<Nz; k++)
	for (j=0; j<Ny; j++)
	  for (i=0; i<Nx; i++)
	    Height.Value(h, k+1, j, i) = Height.Value(h, k, j, i)
	      - r / g * Temperature(h, k, j, i)
	      * log(Pressure(h, k+1, j, i) / Pressure(h, k, j, i));

  }


  //! Computes the altitudes at middle levels from altitudes of interfaces and
  //! pressure and temperature fields.
  /*!
    Level heights are computed according to:
    Z_k = Z_{k-1/2} + a_k * r / g * T_k
    where Z is the altitude, T the temperature, k the level index,
    r the molar gas constant for dry air, g the standard gravity and:
    a_k = 1 - P_{k+1/2} / (P_{k-1/2} - P_{k+1/2}) * log(P_{k-1/2} / P_{k+1/2})
    \param Pressure pressure (Pa).
    \param Temperature temperature (Pa).
    \param InterfHeight (output) altitudes of interfaces (m).
    \param MiddleHeight (output) altitudes of middle points (m).
    \param g (optional) standard gravity. Default: 9.80665.
    \param r (optional) molar gas constant for dry air. Default: 287.0.
    \note Temperature is provided at middle points (not interfaces).
    Pressure is defined at interfaces (including ground-level).
  */
  template<class TP, class TT, class T, class TG>
  void ComputeMiddleHeight(Data<TP, 4, TG>& Pressure,
			   Data<TT, 4, TG>& Temperature,
			   Grid<T>& InterfHeight, Grid<T>& MiddleHeight,
			   T g, T r)
  {

    int h, i, j, k;

    int Nx = MiddleHeight.GetLength(3);
    int Ny = MiddleHeight.GetLength(2);
    int Nz = MiddleHeight.GetLength(1);
    int Nt = MiddleHeight.GetLength(0);

    for (h=0; h<Nt; h++)
      for (k=0; k<Nz; k++)
	for (j=0; j<Ny; j++)
	  for (i=0; i<Nx; i++)
	    MiddleHeight.Value(h, k, j, i) = InterfHeight.Value(h, k, j, i)
	      + r / g * Temperature(h, k, j, i)
	      * (1. - Pressure(h, k+1, j, i)
		 / (Pressure(h, k, j, i) - Pressure(h, k+1, j, i))
		 * log(Pressure(h, k, j, i) / Pressure(h, k+1, j, i)));

  }


  //! Computes the virtual temperature.
  /*!
    The virtual temperature is computed according to: T_v = (1 + c * q) * T
    where T_s is the virtual temperature, q the specific humidity,
    T the temperature and c a coefficient (0.608, usually).
    \param Temperature temperature (K).
    \param SpecificHumidity specific humidity (kg/kg).
    \param VirtualTemperature (output) virtual temperature (K).
    \param c (optional) coefficient. Default: 0.608.
    \note Temperature and VirtualTemperature may be the same object.
  */
  template <class TT, class TH, class T, class TG>
  void ComputeVirtualTemperature(Data<TT, 4, TG>& Temperature,
				 Data<TH, 4, TG>& SpecificHumidity,
				 Data<T, 4, TG>& VirtualTemperature, T c)
  {

    int h, i, j, k;

    int Nx = VirtualTemperature.GetLength(3);
    int Ny = VirtualTemperature.GetLength(2);
    int Nz = VirtualTemperature.GetLength(1);
    int Nt = VirtualTemperature.GetLength(0);

    for (h=0; h<Nt; h++)
      for (k=0; k<Nz; k++)
	for (j=0; j<Ny; j++)
	  for (i=0; i<Nx; i++)
	    VirtualTemperature(h, k, j, i) = Temperature(h, k, j, i)
	      * (1 + c * SpecificHumidity(h, k, j, i));

  }


}  // namespace AtmoData.

#define ATMODATA_FILE_METEOROLOGY_CXX
#endif
