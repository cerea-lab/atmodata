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


#ifndef ATMODATA_FILE_DEPOSITION_HXX

namespace AtmoData
{

  template <class T>
  T ComputeWesely(T surface_temperature, T solar_radiation,
		  string species, T reactivity, T diffusivity, T Henry,
		  T Ri, T Rlu, T Rac, T RgsS, T RgsO, T RclS, T RclO,
		  T limit = 1.e10, T D_H2O = 0.25);

}  // namespace AtmoData.

#define ATMODATA_FILE_DEPOSITION_HXX
#endif
