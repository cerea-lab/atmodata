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


#ifndef ATMODATA_FILE_EMISSIONS_CXX

#include "Emissions.hxx"

namespace AtmoData
{


  //! Computes biogenic emissions.
  /*!
    Computes biogenic emissions according to Simpson et al. (1999).
    \param Temperature soil or leaf temperature.
    \param PAR Photosynthetically active radiation.
    \param LUC Land use coverage.
    \param EF_isoprene isoprene emission factors.
    \param EF_terpenes terpenes emission factors.
    \param EF_NO NO emission factors.
    \param Isoprene isoprene emissions.
    \param Terpenes terpenes emissions.
    \param NO NO emissions.
  */
  template <class TTemp, class TP, class TL, class TD, class TEFI, class TEFT,
	    class TEFN, class TI, class TT, class TN, class TG>
  void ComputeBiogenicEmissions(Data<TTemp, 3, TG>& Temperature, Data<TP, 3, TG>& PAR,
				Data<TL, 3, TG>& LUC, Data<TD, 1, TG>& Density,
				Data<TEFI, 1, TG>& EF_isoprene,
				Data<TEFT, 1, TG>& EF_terpenes,
				Data<TEFN, 1, TG>& EF_NO,
				Data<TI, 3, TG>& Isoprene, Data<TT, 3, TG>& Terpenes, 
				Data<TN, 3, TG>& NO)
  {

    int h, i, j, k;

    int Nx = Isoprene.GetLength(2);
    int Ny = Isoprene.GetLength(1);
    int Nt = Isoprene.GetLength(0);
    int Nc = LUC.GetLength(0);

    TTemp T, Ts_NO;
    TP par;
    double c_l, c_t;

    // For environmental correction factors.
    const double alpha(0.0027), c_l1(1.066);
    const double Ts(303), Tm(314), R(8.314), c_t1(95000), c_t2(230000);
    const double alpha2(alpha * alpha), ratio1(c_t1 / (R * Ts)),
      ratio2(c_t2 / (R * Ts));

    Isoprene.SetZero();
    Terpenes.SetZero();
    NO.SetZero();

    for (h=0; h<Nt; h++)
      for (i=0; i<Nx; i++)
	for (j=0; j<Ny; j++)
	  for (k=0; k<Nc; k++)
	    {
	      // Temperature.
	      T = Temperature(h, j, i);
	      // Photosynthetically active radiation.
	      par = PAR(h, j, i);
	      // Light dependence.
	      c_l = alpha * c_l1 * par
		/ sqrt(1 + alpha2 * par*par);
	      // Temperature dependence.
	      c_t = exp(ratio1 * (T-Ts) / T)
		/ (1. + exp(ratio2 * (T-Tm) / T));
	      // Emission.
	      Isoprene(h, j, i) += Density(k) * EF_isoprene(k)
		* LUC(k, j, i) * c_l * c_t;
	    
	      // Emission.
	      Terpenes(h, j, i) += Density(k) * EF_terpenes(k)
		* LUC(k, j, i) * exp(0.09*(T-Ts));
	    
	      // Emission.
	      if (EF_NO(k)==0.9)
		Ts_NO = .67 * (T-273.15) + 8.8;
	      else
		Ts_NO = .84 * (T-273.15) + 3.6;
	      if (Ts_NO>35)
		Ts_NO = 35;
	      if (Ts_NO<15)
		Ts_NO = 15;
	      NO(h, j, i) += 1.e-3 * EF_NO(k) * LUC(k, j, i) * exp(0.071*(Ts_NO));
	    }

  }


}  // namespace AtmoData.

#define ATMODATA_FILE_EMISSIONS_CXX
#endif
