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
    \param Richardson Richardson number.
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
    \param SurfaceRichardson surface Richardson number.
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
    \param PotentialTemperature potential temperature.
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
    \param Module module.
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
    \param Module surface module.
  */
  template<class TU, class TV, class T, class TG>
  void ComputeModule(Data<TU, 4, TG>& U, Data<TV, 4, TG>& V,
		     Data<T, 3, TG>& Module)
  {

    int h, i, j, k;

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


}  // namespace AtmoData.

#define ATMODATA_FILE_METEOROLOGY_CXX
#endif
