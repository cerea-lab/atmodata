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


#ifndef ATMODATA_FILE_ERRORS_CXX

#include "Errors.hxx"

namespace AtmoData
{

  /*******
   * NGE *
   *******/

  //! Computes the normalized gross error between two data sets.
  /*!
    Data to be compared with reference data is interpolated on
    reference-data grid.
    \param data_ref reference data.
    \param data_comp data to be compared with 'data_ref'.
    \param test boolean function with two parameters; the i-th
    elements of 'data_ref' and 'data_comp' are taken into account
    if 'test(data_ref(i), data_comp(i))' is true.
    \return The normalized gross error.
  */
  template<class T_ref, int N, class TG_ref,
	   class T_comp, class TG_comp>
  T_ref NGE_interpolation(Data<T_ref, N, TG_ref> data_ref,
			  Data<T_comp, N, TG_comp>& data_comp,
			  Function_Base<T_ref, bool>& test)
  {
    Data<T_ref, N, TG_ref> data_comp_interp(data_ref);
    LinearInterpolationGeneral(data_comp, data_comp_interp);

    return NGE(data_ref, data_comp_interp, test);
  }


  //! Computes the normalized gross error between two data sets.
  /*!
    \param data_ref reference data.
    \param data_comp data to be compared with 'data_ref'.
    \param test boolean function with two parameters; the i-th
    elements of 'data_ref' and 'data_comp' are taken into account
    if 'test(data_ref(i), data_comp(i))' is true.
    \return The normalized gross error.
  */
  template<class T_ref, int N, class TG_ref,
	   class T_comp, class TG_comp>
  T_ref NGE(Data<T_ref, N, TG_ref> data_ref,
	    Data<T_comp, N, TG_comp>& data_comp,
	    Function_Base<T_ref, bool>& test)
  {
    T_ref nge;

    T* data_ref_arr = data_ref.GetData();
    T0* data_comp_arr = data_comp.GetData();
    int NbElements = data_ref.GetNbElements();

#ifdef DEBUG_ATMODATA_DIMENSION

    if (NbElements!=data_comp.GetNbElements())
      throw WrongDim("AtmoData::NGE(Data<T_ref, " + to_str(N) +
		     ">&, Data<T_comp, " + to_str(N) +
		     ">&, Function_Base<T_ref, bool>&)",
		     "Data sizes differ.");

#endif
    
    int nb_elt = 0;
    nge = T_ref(0);
    for (int i=0; i<NbElements; i++)
      if (test(data_arr_ref(i), data_arr_comp(i)))
	{
	  nb_elt++;
	  nge += abs( (data_ref_arr[i] - data_comp_arr[i])
		      / data_ref_arr[i] );
	}
    nge = nge / T_ref(nb_elt);

    return nge;
  }


  /********
   * BIAS *
   ********/

  //! Computes the bias between two data sets.
  /*!
    Data to be compared with reference data is interpolated on
    reference-data grid.
    \param data_ref reference data.
    \param data_comp data to be compared with 'data_ref'.
    \param test boolean function with two parameters; the i-th
    elements of 'data_ref' and 'data_comp' are taken into account
    if 'test(data_ref(i), data_comp(i))' is true.
    \return The bias.
  */
  template<class T_ref, int N, class TG_ref,
	   class T_comp, class TG_comp>
  T_ref Bias_interpolation(Data<T_ref, N, TG_ref> data_ref,
			   Data<T_comp, N, TG_comp>& data_comp,
			   Function_Base<T_ref, bool>& test)
  {
    Data<T_ref, N, TG_ref> data_comp_interp(data_ref);
    LinearInterpolationGeneral(data_comp, data_comp_interp);

    return Bias(data_ref, data_comp_interp, test);
  }


  //! Computes the bias between two data sets.
  /*!
    \param data_ref reference data.
    \param data_comp data to be compared with 'data_ref'.
    \param test boolean function with two parameters; the i-th
    elements of 'data_ref' and 'data_comp' are taken into account
    if 'test(data_ref(i), data_comp(i))' is true.
    \return The bias.
  */
  template<class T_ref, int N, class TG_ref,
	   class T_comp, class TG_comp>
  T_ref Bias(Data<T_ref, N, TG_ref> data_ref,
	     Data<T_comp, N, TG_comp>& data_comp,
	     Function_Base<T_ref, bool>& test)
  {
    T_ref bias;

    T* data_ref_arr = data_ref.GetData();
    T0* data_comp_arr = data_comp.GetData();
    int NbElements = data_ref.GetNbElements();

#ifdef DEBUG_ATMODATA_DIMENSION

    if (NbElements!=data_comp.GetNbElements())
      throw WrongDim("AtmoData::Bias(Data<T_ref, " + to_str(N) +
		     ">&, Data<T_comp, " + to_str(N) +
		     ">&, Function_Base<T_ref, bool>&)",
		     "Data sizes differ.");

#endif
    
    int nb_elt = 0;
    bias = T_ref(0);
    for (int i=0; i<NbElements; i++)
      if (test(data_arr_ref(i), data_arr_comp(i)))
	{
	  nb_elt++;
	  bias += data_ref_arr[i] - data_comp_arr[i];
	}
    bias = bias / T_ref(nb_elt);

    return bias;
  }


  /********
   * RMS *
   ********/

  //! Computes the root mean square between two data sets.
  /*!
    Data to be compared with reference data is interpolated on
    reference-data grid.
    \param data_ref reference data.
    \param data_comp data to be compared with 'data_ref'.
    \param test boolean function with two parameters; the i-th
    elements of 'data_ref' and 'data_comp' are taken into account
    if 'test(data_ref(i), data_comp(i))' is true.
    \return The root mean square.
  */
  template<class T_ref, int N, class TG_ref,
	   class T_comp, class TG_comp>
  T_ref RMS_interpolation(Data<T_ref, N, TG_ref> data_ref,
			  Data<T_comp, N, TG_comp>& data_comp,
			  Function_Base<T_ref, bool>& test)
  {
    Data<T_ref, N, TG_ref> data_comp_interp(data_ref);
    LinearInterpolationGeneral(data_comp, data_comp_interp);

    return RMS(data_ref, data_comp_interp, test);
  }


  //! Computes the root mean square between two data sets.
  /*!
    \param data_ref reference data.
    \param data_comp data to be compared with 'data_ref'.
    \param test boolean function with two parameters; the i-th
    elements of 'data_ref' and 'data_comp' are taken into account
    if 'test(data_ref(i), data_comp(i))' is true.
    \return The root mean square.
  */
  template<class T_ref, int N, class TG_ref,
	   class T_comp, class TG_comp>
  T_ref RMS(Data<T_ref, N, TG_ref> data_ref,
	    Data<T_comp, N, TG_comp>& data_comp,
	    Function_Base<T_ref, bool>& test)
  {
    T_ref rms;

    T* data_ref_arr = data_ref.GetData();
    T0* data_comp_arr = data_comp.GetData();
    int NbElements = data_ref.GetNbElements();

#ifdef DEBUG_ATMODATA_DIMENSION

    if (NbElements!=data_comp.GetNbElements())
      throw WrongDim("AtmoData::RMS(Data<T_ref, " + to_str(N) +
		     ">&, Data<T_comp, " + to_str(N) +
		     ">&, Function_Base<T_ref, bool>&)",
		     "Data sizes differ.");

#endif
    
    int nb_elt = 0;
    bias = T_ref(0);
    for (int i=0; i<NbElements; i++)
      if (test(data_arr_ref(i), data_arr_comp(i)))
	{
	  nb_elt++;
	  rms += (data_arr[i] - data0_arr[i]) * (data_arr[i] - data0_arr[i]);
	}
    rms = rms / T_ref(nb_elt);

    return rms;
  }


}  // namespace AtmoData.


#define ATMODATA_FILE_ERRORS_CXX
#endif
