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


#ifndef ATMODATA_FILE_EMISSIONS_HXX

namespace AtmoData
{

  template <class TL, class TD, class TEFI, class TEFT,
	    class TEFN, class TI, class TT, class TN, class TG>
  void ComputeBiogenicRates(Data<TL, 3, TG>& LUC, Data<TD, 3, TG>& Density,
			    Data<TEFI, 1, TG>& EF_isoprene,
			    Data<TEFT, 1, TG>& EF_terpenes,
			    Data<TEFN, 1, TG>& EF_NO,
			    Data<TL, 3, TG>& Isoprene, Data<TL, 3, TG>& Terpenes, 
			    Data<TL, 3, TG>& NO);

  template <class TTemp, class TP, class TL, class TD, class TEFI, class TEFT,
	    class TEFN, class TI, class TT, class TN, class TG>
  void ComputeBiogenicEmissions(Data<TTemp, 3, TG>& Temperature, Data<TP, 3, TG>& PAR,
				Data<TL, 3, TG>& LUC, Data<TD, 3, TG>& Density,
				Data<TEFI, 1, TG>& EF_isoprene,
				Data<TEFT, 1, TG>& EF_terpenes,
				Data<TEFN, 1, TG>& EF_NO,
				Data<TL, 3, TG>& Isoprene, Data<TL, 3, TG>& Terpenes, 
				Data<TL, 3, TG>& NO);

}  // namespace AtmoData.

#define ATMODATA_FILE_EMISSIONS_HXX
#endif
