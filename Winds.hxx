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


#ifndef ATMODATA_FILE_WINDS_HXX

namespace AtmoData
{

  template<class TU, class TV, class TW, class TG>
  void GetVerticalWind(Data<TU, 4, TG>& U, Data<TV, 4, TG>& V, Data<TW, 4, TG>& W);

}  // namespace AtmoData.

#define ATMODATA_FILE_WINDS_HXX
#endif
