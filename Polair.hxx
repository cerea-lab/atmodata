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


#ifndef ATMODATA_FILE_POLAIR_HXX

namespace AtmoData
{

  template<class T, class TG>
  void TransformZonalWind(Data<T, 4, TG>& ZonalWind);

  template<class T, class TG>
  void TransformMeridionalWind(Data<T, 4, TG>& MeridionalWind);

}  // namespace AtmoData.

#define ATMODATA_FILE_POLAIR_HXX
#endif
