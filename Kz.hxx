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


#ifndef ATMODATA_FILE_KZ_HXX

namespace AtmoData
{

  //! Returns the minimum of two elements.
  /*! \param x first number.
    \param y second number.
    \return The minimum of {x, y}.
  */
  template <class T, class T0>
  inline T min(T x, T0 y)
  {

    return x<y?x:y;

  }

  //! Returns the maximum of two elements.
  /*! \param x first number.
    \param y second number.
    \return The maximum of {x, y}.
  */
  template <class T, class T0>
  inline T max(T x, T0 y)
  {

    return x>y?x:y;

  }

  template<class TU, class TV, class TW,
	   class TTp, class T, class TG>
  void ComputeLouisKz(Data<TU, 4, TG>& U, Data<TV, 4, TG>& V, Data<TW, 4, TG>& W,
		      Data<TTp, 4, TG>& Tp, Data<T, 4, TG>& Kz,
		      T B = T(5), T C = T(5), T D = T(5),
		      T Ka = T(0.4), T z0 = T(1), T L0 = T(100));

}  // namespace AtmoData.

#define ATMODATA_FILE_KZ_HXX
#endif
