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
  void LouisKz(Data<TU, 4, TG>& U, Data<TV, 4, TG>& V, Data<TW, 4, TG>& W,
	       Data<TTp, 4, TG>& Tp, Data<T, 4, TG>& Kz,
	       T B = T(5), T C = T(5), T D = T(5),
	       T Ka = T(0.4), T z0 = T(1), T L0 = T(100));

}  // namespace AtmoData.

#define ATMODATA_FILE_KZ_HXX
#endif
