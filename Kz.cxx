// Copyright (C) 2003 Vivien Mallet
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

  //! Computes vertical diffusion according to the Louis formula.
  /*!
    \param U zonal wind.
    \param V meridional wind.
    \param W vertical wind.
    \param Tp potential temperature.
    \param Kz vertical diffusion coefficients.
    \param B parameter.
    \param C parameter.
    \param D parameter.
    \param Ka Von Karman constant.
    \param z0 scale parameter.
    \param L0 scale parameter.
  */
  template<class TU, class TV, class TW,
	   class TTp, class T, class TG>
  void LouisKz(Data<TU, 4, TG>& U, Data<TV, 4, TG>& V, Data<TW, 4, TG>& W,
	       Data<TTp, 4, TG>& Tp, Data<T, 4, TG>& Kz,
	       T B, T C, T D, T Ka, T z0, T L0)
  {

    int h, i, j, k;

    int Nx = W.GetLength(3);
    int Ny = W.GetLength(2);
    int Nz = W.GetLength(1) - 1;
    int Nt = min( min(U.GetLength(0), V.GetLength(0)), W.GetLength(0) );

    T l, R, F, L;
    T dWind_dz, derivative;
    T g = 9.81;

    Grid<T>& Levels = W[1];
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

	      // dW/dz.
	      
	      derivative = ( W(h, k+1, j, i) - W(h, k-1, j, i) )
		/ ( Levels(k+1) - Levels(k) ) * 0.5;
	    
	      derivative = derivative * derivative;
	      dWind_dz = derivative;

	      // dU/dz.

	      derivative = ( U(h, k, j, i+1) + U(h, k, j, i)
			     - U(h, k-1, j, i+1) - U(h, k-1, j, i) )
		/ ( Nodes(k) - Nodes(k-1) ) * 0.5;

	      derivative = derivative * derivative;
	      dWind_dz += derivative;

	      // dV/dz.

	      derivative = ( V(h, k, j+1, i) + V(h, k, j, i)
			     - V(h, k-1, j+1, i) - V(h, k-1, j, i) )
		/ ( Nodes(k) - Nodes(k-1) ) * 0.5;
	      
	      derivative = derivative * derivative;
	      dWind_dz += derivative;

	      dWind_dz = sqrt(dWind_dz);

	      /* dWind / dz */
	      /**************/

	      /***********/
	      /* F(R, z) */

	      derivative = log( (Tp(h, k, j, i)) / (Tp(h, k-1, j, i)) )
		/ ( Nodes(k) - Nodes(k-1) );

	      R = g * derivative / ( dWind_dz * dWind_dz );

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
