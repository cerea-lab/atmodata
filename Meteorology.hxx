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


#ifndef ATMODATA_FILE_METEOROLOGY_HXX

namespace AtmoData
{

  template<class TU, class TV, class TTp, class T, class TG>
  void ComputeRichardson(Data<TU, 4, TG>& ZonalWind,
			 Data<TV, 4, TG>& MeridionalWind,
			 Data<TTp, 4, TG>& PotentialTemperature,
			 Data<T, 4, TG>& Richardson,
			 T wind_threshold = 0.001);

  template<class TU, class TV, class TTp, class T, class TG>
  void ComputeRichardson(Data<TU, 4, TG>& ZonalWind,
			 Data<TV, 4, TG>& MeridionalWind,
			 Data<TTp, 4, TG>& PotentialTemperature,
			 Data<T, 3, TG>& SurfaceRichardson,
			 T wind_threshold = 0.001);

  template<class TU, class TTp, class T, class TG>
  void ComputeRichardson(Data<TU, 3, TG>& WindModule,
			 Data<TTp, 3, TG>& SurfacePotentialTemperature,
			 Data<TTp, 4, TG>& PotentialTemperature,
			 Data<T, 3, TG>& SurfaceRichardson,
			 T wind_threshold = 0.001);

  template<class TR, class TU, class TTp, class T, class TG>
  void ComputeRichardson(Data<TR, 2, TG>& Roughness,
			 Data<TU, 3, TG>& WindModule,
			 Data<TTp, 3, TG>& SurfacePotentialTemperature,
			 Data<TTp, 4, TG>& PotentialTemperature,
			 Data<T, 3, TG>& SurfaceRichardson,
			 T wind_threshold = 0.001);

  template<class TT, class TP, class T, class TG>
  void ComputePotentialTemperature(Data<TT, 4, TG>& Temperature,
				   Data<TP, 4, TG>& Pressure,
				   Data<T, 4, TG>& PotentialTemperature,
				   T P0 = 101325., T cp = 1005.,T r = 287.0);

  template<class TT, class TP, class T, class TG>
  void ComputePotentialTemperature(Data<TT, 3, TG>& Temperature,
				   Data<TP, 3, TG>& Pressure,
				   Data<T, 3, TG>& PotentialTemperature,
				   T P0 = 101325., T cp = 1005.,T r = 287.0);

  template<class TT, class TP, class T, class TG>
  void ComputeSaturationHumidity(Data<TT, 3, TG>& Temperature,
				 Data<TP, 3, TG>& Pressure,
				 Data<T, 3, TG>& SaturationHumidity);

  template<class TT, class TP, class T, class TG>
  void ComputeSaturationHumidity(Data<TT, 4, TG>& Temperature,
				 Data<TP, 4, TG>& Pressure,
				 Data<T, 4, TG>& SaturationHumidity);

  template<class TS, class TT, class TP, class T, class TG>
  void ComputeRelativeHumidity(Data<TS, 4, TG>& SpecificHumidity,
			       Data<TT, 4, TG>& Temperature,
			       Data<TP, 4, TG>& Pressure,
			       Data<T, 4, TG>& RelativeHumidity);

  template<class TH, class TS, class TW, class TL, class T, class TG>
  void ComputeSurfaceHumidity_diag(Data<TH, 4, TG>& Humidity,
				   Data<TS, 3, TG>& SaturationHumidity,
				   Data<TW, 3, TG>& SoilWater,
				   Data<TL, 3, TG>& LUC, int sea_index,
				   Data<T, 3, TG>& SurfaceHumidity,
				   T veg = 1.0, T theta_cap = 0.323);

  template<class TS, class TP, class T, class TG>
  void ComputeCriticalRelativeHumidity(Data<TS, 3, TG>& SurfacePressure,
				       Data<TP, 4, TG>& Pressure,
				       Data<T, 4, TG>& CriticalRelativeHumidity,
				       T coeff0 = 2., T coeff1 = sqrt(3.));

  template<class TS, class TP, class T, class TG>
  void ComputeCriticalRelativeHumidity_extended(Data<TS, 3, TG>& SurfacePressure,
						Data<TP, 4, TG>& Pressure,
						Data<T, 4, TG>& CriticalRelativeHumidity,
						T coeff0 = 1.1,
						T coeff1 = sqrt(1.3),
						T a0 = 0., T a1 = 1.1);

  template<class TB, class TS, class TP, class T, class TG>
  void ComputeCriticalRelativeHumidity(Data<TB, 3, TG>& BoundaryLayerHeight,
				       Data<TS, 3, TG>& SurfacePressure,
				       Data<TP, 4, TG>& Pressure,
				       Data<T, 4, TG>& CriticalRelativeHumidity,
				       T coeff0 = 2., T coeff1 = sqrt(3.),
				       T BL_CRH = 0.98);

  template<class TP, class T, class TG>
  void ComputeCriticalRelativeHumidity(Data<TP, 4, TG>& Pressure,
				       Data<T, 4, TG>& CriticalRelativeHumidity,
				       T CRH_0 = 0.75, T CRH_1 = 0.95,
				       T CRH_2 = 0.95,
				       T P_0 = 70000., T P_1 = 40000.);

  template<class TR, class TC, class T, class TG>
  void ComputeCloudFraction(Data<TR, 4, TG>& RelativeHumidity,
			    Data<TC, 4, TG>& CriticalRelativeHumidity,
			    Data<T, 4, TG>& CloudFraction);

  template<class TP, class TR, class TC, class T, class TG>
  void ComputeCloudFraction(Data<TP, 3, TG>& BoundaryLayerHeight,
			    Data<TR, 4, TG>& RelativeHumidity,
			    Data<TC, 4, TG>& CriticalRelativeHumidity,
			    Data<T, 4, TG>& CloudFraction);

  template<class TU, class TV, class T, class TG>
  void ComputeModule(Data<TU, 4, TG>& U, Data<TV, 4, TG>& V,
		     Data<T, 4, TG>& Module);

  template<class TU, class TV, class T, class TG>
  void ComputeModule(Data<TU, 4, TG>& U, Data<TV, 4, TG>& V,
		     Data<T, 3, TG>& Module);

  template<class TU, class TV, class T, class TG>
  void ComputeModule(Data<TU, 3, TG>& U, Data<TV, 3, TG>& V,
		     Data<T, 3, TG>& Module);

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
			 T P_0 = 80000., T P_1 = 45000.);

  template <class TT, class TP, class TH,
	    class T, class TG>
  void ComputeCloudHeight(Data<TT, 4, TG>& Temperature,
			  Data<TP, 4, TG>& Pressure,
			  Data<TH, 4, TG>& Humidity,
			  T (CriticalRelativeHumidity)(const T&, const T&,
						       const T&),
			  Data<T, 3, TG>& CloudHeight);

  template <class TT, class TP, class TH,
	    class TCRH, class T, class TG>
  void ComputeCloudHeight(Data<TT, 4, TG>& Temperature,
			  Data<TP, 4, TG>& Pressure,
			  Data<TH, 4, TG>& Humidity,
			  Data<TCRH, 4, TG>& CriticalRelativeHumidity,
			  Data<T, 3, TG>& CloudHeight);

  template <class T, class TG>
  void ComputeCloudHeight(Data<int, 4>& LowIndices,
			  Data<int, 4>& MediumIndices,
			  Data<int, 4>& HighIndices,
			  Grid<TG>& GridZ_interf,
			  Data<T, 3, TG>& CloudHeight);

  template <class Ta, class Tb, class TSP,
	    class T, class TG>
  void ComputePressure(Data<Ta, 1, TG>& alpha, Data<Tb, 1, TG>& beta,
		       Data<TSP, 3, TG>& SurfacePressure,
		       Data<T, 4, TG>& Pressure, T P0 = 101325.);

  template<class TPS, class TP, class TT, class T, class TG>
  void ComputeHeight(Data<TPS, 3, TG>& SurfacePressure,
		     Data<TP, 4, TG>& Pressure,
		     Data<TT, 4, TG>& Temperature,
		     Grid<T>& Height, T g = 9.80665, T r = 287.0);

  template<class TP, class TT, class T, class TG>
  void ComputeInterfHeight(Data<TP, 4, TG>& Pressure,
			   Data<TT, 4, TG>& Temperature,
			   Grid<T>& Height, bool ground_set = false,
			   T g = 9.80665, T r = 287.0);

  template<class TP, class TT, class T, class TG>
  void ComputeMiddleHeight(Data<TP, 4, TG>& Pressure,
			   Data<TT, 4, TG>& Temperature,
			   Grid<T>& InterfHeight, Grid<T>& MiddleHeight,
			   T g = 9.80665, T r = 287.0);

  template <class TT, class TH, class T, class TG>
  void ComputeVirtualTemperature(Data<TT, 4, TG>& Temperature,
				 Data<TH, 4, TG>& SpecificHumidity,
				 Data<T, 4, TG>& VirtualTemperature,
				 T c = 0.608);

}  // namespace AtmoData.

#define ATMODATA_FILE_METEOROLOGY_HXX
#endif
