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


#ifndef ATMODATA_FILE_METEOROLOGY_CXX

#include "Meteorology.hxx"

namespace AtmoData
{


  //! Computes the Richardson number.
  /*!
    Winds may be provided in two ways. The first option is to provide winds on interfaces
    (along x for the zonal wind, along y for the meridional wind). The second option is
    simply to provide winds at nodes (i.e. where the potential temperature and the Richardson
    number are defined).
    \param ZonalWind zonal wind.
    \param MeridionalWind meridional wind.
    \param PotentialTemperature potential temperature.
    \param Richardson (output) Richardson number.
    \param wind_threshold (optional) minimum of the wind module. Default: 0.1.
  */
  template<class TU, class TV, class TTp, class T, class TG>
  void ComputeRichardson(Data<TU, 4, TG>& ZonalWind, Data<TV, 4, TG>& MeridionalWind,
			 Data<TTp, 4, TG>& PotentialTemperature, Data<T, 4, TG>& Richardson,
			 T wind_threshold)
  {

    int h, i, j, k;

    int Nx = Richardson.GetLength(3);
    int Ny = Richardson.GetLength(2);
    int Nz = Richardson.GetLength(1);
    int Nt = Richardson.GetLength(0);

    Grid<TG>& Levels = Richardson[1];

    const T g(9.81);
    T u, v, wind;
    int level;

    if ( (ZonalWind.GetLength(3) == Nx + 1) && (MeridionalWind.GetLength(2) == Ny + 1) )
      for (h=0; h<Nt; h++)
	for (k=0; k<Nz; k++)
	  for (j=0; j<Ny; j++)
	    for (i=0; i<Nx; i++)
	      {
		
		level = k==0 ? 1 : k;
		
		u = 0.5 * (ZonalWind(h, k, j, i+1) + ZonalWind(h, k, j, i));
		v = 0.5 * (MeridionalWind(h, k, j+1, i) + MeridionalWind(h, k, j, i));
		wind = max(sqrt(u*u + v*v), wind_threshold);
		SurfaceRichardson(h, k, j, i) = g * (PotentialTemperature(h, level, j, i)
						     - PotentialTemperature(h, level - 1, j, i))
		  * (Levels.Value(h, level, j, i) - Levels.Value(h, level - 1, j, i))
		  / (wind*wind * PotentialTemperature(h, level - 1, j, i));

	      }
    else
      for (h=0; h<Nt; h++)
	for (k=0; k<Nz; k++)
	  for (j=0; j<Ny; j++)
	    for (i=0; i<Nx; i++)
	      {
		
		level = k==0 ? 1 : k;
		
		u = ZonalWind(h, k, j, i);
		v = MeridionalWind(h, k, j, i);
		wind = max(sqrt(u*u + v*v), wind_threshold);
		SurfaceRichardson(h, k, j, i) = g * (PotentialTemperature(h, level, j, i)
						     - PotentialTemperature(h, level - 1, j, i))
		  * (Levels.Value(h, level, j, i) - Levels.Value(h, level - 1, j, i))
		  / (wind*wind * PotentialTemperature(h, level - 1, j, i));

	      }

  }


  //! Computes the surface Richardson number.
  /*!
    Winds may be provided in two ways. The first option is to provide winds on interfaces
    (along x for the zonal wind, along y for the meridional wind). The second option is
    simply to provide winds at nodes (i.e. where the potential temperature and the Richardson
    number are defined).
    \param ZonalWind zonal wind.
    \param MeridionalWind meridional wind.
    \param PotentialTemperature potential temperature.
    \param SurfaceRichardson (output) surface Richardson number.
    \param wind_threshold (optional) minimum of the wind module. Default: 0.1.
  */
  template<class TU, class TV, class TTp, class T, class TG>
  void ComputeRichardson(Data<TU, 4, TG>& ZonalWind, Data<TV, 4, TG>& MeridionalWind,
			 Data<TTp, 4, TG>& PotentialTemperature,
			 Data<T, 3, TG>& SurfaceRichardson,
			 T wind_threshold)
  {

    int h, i, j;

    int Nx = SurfaceRichardson.GetLength(2);
    int Ny = SurfaceRichardson.GetLength(1);
    int Nt = SurfaceRichardson.GetLength(0);

    Grid<TG>& Levels = PotentialTemperature[1];

    const T g(9.81);
    T u, v, wind;

    if ( (ZonalWind.GetLength(3) == Nx + 1) && (MeridionalWind.GetLength(2) == Ny + 1) )
      for (h=0; h<Nt; h++)
	for (j=0; j<Ny; j++)
	  for (i=0; i<Nx; i++)
	    {
		
	      u = 0.5 * (ZonalWind(h, 0, j, i+1) + ZonalWind(h, 0, j, i));
	      v = 0.5 * (MeridionalWind(h, 0, j+1, i) + MeridionalWind(h, 0, j, i));
	      wind = max(sqrt(u*u + v*v), wind_threshold);
	      SurfaceRichardson(h, j, i) = g * (PotentialTemperature(h, 1, j, i)
						- PotentialTemperature(h, 0, j, i))
		* (Levels.Value(h, 1, j, i) - Levels.Value(h, 0, j, i))
		/ (wind*wind * PotentialTemperature(h, 0, j, i));

	    }
    else
      for (h=0; h<Nt; h++)
	for (j=0; j<Ny; j++)
	  for (i=0; i<Nx; i++)
	    {
		
	      u = ZonalWind(h, 0, j, i);
	      v = MeridionalWind(h, 0, j, i);
	      wind = max(sqrt(u*u + v*v), wind_threshold);
	      SurfaceRichardson(h, j, i) = g * (PotentialTemperature(h, 1, j, i)
						- PotentialTemperature(h, 0, j, i))
		* (Levels.Value(h, 1, j, i) - Levels.Value(h, 0, j, i))
		/ (wind*wind * PotentialTemperature(h, 0, j, i));

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
    \param r (optional) molar gas constant for air. Default: 287.0 J.kg^{-1}.K^{-1}.
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


  //! Computes the module of a 2D-vectors field.
  /*!
    This function was initially dedicated to winds. In this case, zonal winds and meridional
    winds are provided and the module of the wind is computed (assuming that the vertical
    wind is zero). Winds may be provided in two ways. The first option is to provide winds
    on interfaces (along x for the zonal wind, along y for the meridional wind). 
    The second option is simply to provide winds at nodes (i.e. where the module is defined).
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
    This function was initially dedicated to winds. In this case, zonal winds and meridional
    winds are provided and the module of the wind is computed (assuming that the vertical
    wind is zero). Winds may be provided in two ways. The first option is to provide winds
    on interfaces (along x for the zonal wind, along y for the meridional wind). 
    The second option is simply to provide winds at nodes (i.e. where the module is defined).
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
    This function was initially dedicated to winds. In this case, zonal winds and meridional
    winds are provided and the module of the wind is computed (assuming that the vertical
    wind is zero). Winds may be provided in two ways. The first option is to provide winds
    on interfaces (along x for the zonal wind, along y for the meridional wind). 
    The second option is simply to provide winds at nodes (i.e. where the module is defined).
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


  //! Computes the height of cloud basis.
  /*!
    \param Temperature temperature (K).
    \param Pressure pressure (Pa).
    \param Humidity specific humidity (kg/kg).
    \param CriticalRelativeHumidity function that returns the critical
    relative humidity as function of the altitude, the pressure and reference pressure.
    \param CloudHeight (output) altitudes of cloud basis.
  */
  template <class TT, class TP, class TH,
	    class T, class TG>
  void ComputeCloudHeight(Data<TT, 4, TG>& Temperature, Data<TP, 4, TG>& Pressure,
			  Data<TH, 4, TG>& Humidity,
			  T (CriticalRelativeHumidity)(const T&, const T&, const T&),
			  Data<T, 3, TG>& CloudHeight)
  {

    int h, k, j, i;
    int Nt(CloudHeight.GetLength(0));
    int Nz(Pressure.GetLength(1));
    int Ny(CloudHeight.GetLength(1));
    int Nx(CloudHeight.GetLength(2));    

    // Index "0" and "1" refer to two contiguous levels.
    T rh0, rh1, rhc, dz, delta_z, s;

    T max_height = 2. * Pressure[1].Value(0, Nz-1, 0, 0);
    CloudHeight.Fill(max_height);

    for (h=0; h<Nt; h++)
      for (j=0; j<Ny; j++)
	for (i=0; i<Nx; i++)
	  {

	    // Specific humidity to relative humidity.
	    s = 611. * pow(10., 7.5 * (Temperature(h, Nz-1, j, i) - 273.15)
			   / (Temperature(h, Nz-1, j, i) - 35.85));
	    rh0 = Pressure(h, Nz-1, j, i) * Humidity(h, Nz-1, j, i)
	      / (0.62197 + Humidity(h, Nz-1, j, i)) / s;
	    
	    k = 0;
	    while ( (k<Nz) && (CloudHeight(h, j, i) == max_height) )
	      {

		// Specific humidity to relative humidity.
		s = 611. * pow(10., 7.5 * (Temperature(h, k, j, i) - 273.15)
			       / (Temperature(h, k, j, i) - 35.85));
		rh1 = Pressure(h, k, j, i) * Humidity(h, k, j, i)
		  / (0.62197 + Humidity(h, k, j, i)) / s;

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
  void ComputeHeight(Data<TPS, 3, TG>& SurfacePressure, Data<TP, 4, TG>& Pressure,
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


  //! Computes the altitudes at interfaces from pressure and temperature fields.
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
    false if they have to be set (they are set to zero in this case). Default: false.
    \note Temperature is provided at middle points (not interfaces).
    Pressure and Height are defined at interfaces (including ground-level).
  */
  template<class TP, class TT, class T, class TG>
  void ComputeInterfHeight(Data<TP, 4, TG>& Pressure, Data<TT, 4, TG>& Temperature,
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
  void ComputeMiddleHeight(Data<TP, 4, TG>& Pressure, Data<TT, 4, TG>& Temperature,
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
	      - r / g * Temperature(h, k, j, i)
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
