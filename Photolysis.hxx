// Copyright (C) 2003-2004 Ecole Nationale des Ponts et Chaussees
//     Author: Vivien Mallet
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


#ifndef ATMODATA_FILE_PHOTOLYSIS_HXX

namespace AtmoData
{

  template<class T>
  T ZenithAngle(T lon, T lat, int idate, T ut);

  template <class TH, class TL, class TMC, class THC, class T, class TG>
  void ComputeAttenuation_LWC(Data<TH, 4, TG>& Humidity, Data<TH, 4, TG>& CRH,
			      Data<TL, 4, TG>& LiquidWaterContent,
			      Data<TMC, 3, TG>& MediumCloudiness,
			      Data<THC, 3, TG>& HighCloudiness,
			      int date, Data<T, 4, TG>& Attenuation);

  template <class TL, class TMC, class THC, class T, class TG>
  void ComputeAttenuation_LWC(Data<TL, 4, TG>& LiquidWaterContent,
			      Data<int, 4> LowIndices,
			      Data<int, 4> MediumIndices,
			      Data<int, 4> HighIndices,
			      Data<TMC, 3, TG>& MediumCloudiness,
			      Data<THC, 3, TG>& HighCloudiness,
			      int date, Data<T, 4, TG>& Attenuation);

}  // namespace AtmoData.

#define ATMODATA_FILE_PHOTOLYSIS_HXX
#endif
