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


#ifndef ATMODATA_FILE_KZ_CXX

#include "Kz.hxx"

namespace AtmoData
{

  //! Computes vertical diffusion according to the Louis formula (1979).
  /*!
    \param U zonal wind.
    \param V meridional wind.
    \param Tp potential temperature.
    \param Kz (output) vertical diffusion coefficients at the interfaces.
    \param L0 scale parameter. Default: 100.
    \param B parameter. Default: 5.
    \param C parameter. Default: 5.
    \param D parameter. Default: 5.
    \param z0 scale parameter. Default: 1.
    \param a parameter. Default: 0.115.
    \param b parameter. Default: 0.175.
    \param delta_z0 parameter. Default: 0.01.
    \param Ka Von Karman constant. Default: 0.4.
    \note Kz is given at the interfaces. It is assumed to be 0 at the first
    interface and at the last one.
    \note The critical Richardson number is computed following Pielke (2002)
    and Nordeng (1986): Ric = a * (delta_z / delta_z0)^b.
  */
  template<class TU, class TV, class TTp, class T, class TG>
  void ComputeLouisKz(Data<TU, 4, TG>& U, Data<TV, 4, TG>& V,
		      Data<TTp, 4, TG>& Tp, Data<T, 4, TG>& Kz,
		      T L0, T B, T C, T D, T z0, T a, T b, T delta_z0, T Ka)
  {

    int h, i, j, k;

    int Nx = Kz.GetLength(3);
    int Ny = Kz.GetLength(2);
    int Nz = Kz.GetLength(1) - 1;
    int Nt = min(U.GetLength(0), V.GetLength(0));

    T l, R, F, L;
    T dWind_dz, derivative, dz;
    T g = 9.81;

    Grid<T>& Levels = Kz[1];
    Grid<T>& Nodes = U[1];

    Kz.SetZero();

    for (h=0; h<Nt; h++)
      for (k=1; k<Nz; k++)
	for (j=0; j<Ny; j++)
	  for (i=0; i<Nx; i++)
	    {

	      /*********************/
	      /* l = Ka * ---- ... */
	      
	      L = L0;
	      
	      l = Ka * (Levels(k)+z0) / (1.0 + Ka * (Levels(k)+z0) / L);
	      
	      /* l = Ka * ---- ... */
	      /*********************/

	      /**************/
	      /* dWind / dz */

	      dz = Nodes(k) - Nodes(k-1);

	      // dU/dz.

	      derivative = ( U(h, k, j, i+1) + U(h, k, j, i)
			     - U(h, k-1, j, i+1) - U(h, k-1, j, i) )
		/ dz * 0.5;

	      derivative = derivative * derivative;
	      dWind_dz = derivative;

	      // dV/dz.

	      derivative = ( V(h, k, j+1, i) + V(h, k, j, i)
			     - V(h, k-1, j+1, i) - V(h, k-1, j, i) )
		/ dz * 0.5;
	      
	      derivative = derivative * derivative;
	      dWind_dz += derivative;

	      dWind_dz = sqrt(dWind_dz);

	      /* dWind / dz */
	      /**************/

	      /***********/
	      /* F(R, z) */

              derivative = ( Tp(h, k, j, i) - Tp(h, k-1, j, i) )
                / (dz * Tp(h, k-1, j, i));

	      R = min(g * derivative / ( dWind_dz * dWind_dz ),
		      a * pow(dz / delta_z0, b));

	      if (R>=0)
		F = 1.0 / ( 1.0 + 3.0 * B * R * sqrt(1.0+D*R) );
	      else
		{
		  F = 1.0 + 3.0 * B * C * sqrt(fabs(R)/27.0) *
		    (l * l) / ( (Levels(k) + z0) * (Levels(k) + z0) );
		  F = 1.0 - 3.0 * B * R / F;
		}

	      /* F(R, z) */
	      /***********/

	      Kz(h, k, j, i) = l * l * dWind_dz * F;

	    }

  }


}  // namespace AtmoData.

#define ATMODATA_FILE_KZ_CXX
#endif
