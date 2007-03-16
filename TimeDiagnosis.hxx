// Copyright (C) 2006 CEREA
//     Author: Irène Korsakissok
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
//     http://www.enpc.fr/cerea/atmodata/


#ifndef ATMODATA_FILE_TIMEDIAGNOSIS_HXX


namespace AtmoData
{
  template<class T>
  void ComputeDeclination(Date date, T& declination, T& time_equation);

  template<class T>
  void ComputeDeclination(int idate, T ut, T& declination, T& time_equation);

  template<class T>
  void ComputeSunHour(T lon, T lat, int idate,
		      T& sunrise_hour, T& sunset_hour);

  template<class T>
  T ComputeSunriseHour(T lon, T lat, int idate);

  template<class T>
  T ComputeSunsetHour(T lon, T lat, int idate);

  template<class T>
  bool IsDay(T lon, T lat, int date, T ut);

  template<class T>
  bool IsDay(T lon, T lat, Date date);


}  // namespace AtmoData.


#define ATMODATA_FILE_TIMEDIAGNOSIS_HXX
#endif
