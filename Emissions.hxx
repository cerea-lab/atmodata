// Copyright (C) 2003-2005 CEREA
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


#ifndef ATMODATA_FILE_EMISSIONS_HXX

#include <list>

namespace AtmoData
{

  template <class TL, class TD, class TEFI, class TEFT,
	    class TEFN, class TI, class TT, class TN, class TG>
  void ComputeBiogenicRates(Data<TL, 3, TG>& LUC, Data<TD, 3, TG>& Density,
			    Data<TEFI, 1, TG>& EF_isoprene,
			    Data<TEFT, 1, TG>& EF_terpenes,
			    Data<TEFN, 1, TG>& EF_NO,
			    Data<TL, 3, TG>& Isoprene,
			    Data<TL, 3, TG>& Terpenes, 
			    Data<TL, 3, TG>& NO);

  template <class TTemp, class TP, class TL, class TD, class TEFI, class TEFT,
	    class TEFN, class TI, class TT, class TN, class TG>
  void ComputeBiogenicEmissions(Data<TTemp, 3, TG>& Temperature,
				Data<TP, 3, TG>& PAR,
				Data<TL, 3, TG>& LUC,
				Data<TD, 3, TG>& Density,
				Data<TEFI, 1, TG>& EF_isoprene,
				Data<TEFT, 1, TG>& EF_terpenes,
				Data<TEFN, 1, TG>& EF_NO,
				Data<TL, 3, TG>& Isoprene,
				Data<TL, 3, TG>& Terpenes, 
				Data<TL, 3, TG>& NO);

  //! Stores an EMEP emission associated with a given country.
  template <class T>
  class EmepCountryEmission
  {
  public:
    //! Emission.
    T emission_;
    //! EMEP country number.
    int country_;
  public:
    EmepCountryEmission(T emission, int country);
  };

  //! Provides time zone offset from GMT for a list of countries.
  class TimeZone
  {

  public:
    //! EMEP country number.
    vector<int> countries_;
    //! Time zone offset from GMT.
    vector<int> local_times_;

    TimeZone(int N);
    void Init(string file_name);
    int operator () (int i) const;
  };


}  // namespace AtmoData.

#define ATMODATA_FILE_EMISSIONS_HXX
#endif
