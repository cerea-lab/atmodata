#ifndef ATMODATA_FILE_WINDS_CXX

#include "Winds.hxx"

namespace AtmoData
{

  //! Computes vertical wind-speeds from zonal and meridional wind-speeds.
  /*! Coordinates: W(t, z, y, x).
    \param U zonal wind.
    \param V meridional wind.
    \param W vertical wind.
  */
  template<class TU, class TV, class TW, class TG>
  void GetVerticalWind(Data<TU, 4, TG>& U, Data<TV, 4, TG>& V, Data<TW, 4, TG>& W)
  {
  
    int h, i, j, k;

    const TW pi = 3.14159265358979323846264;

    int Nt = W.GetLength(0);
    int Nz = W.GetLength(1);
    int Ny = W.GetLength(2);
    int Nx = W.GetLength(3);

    TW R_Earth = 6371229.;

    for (h=0; h<Nt; h++)
      for (j=0; j<Ny; j++)
	for (i=0; i<Nx; i++)
	  W(h, 0, j, i) = TW(0);
  
    for (h=0; h<Nt; h++)
      for (k=1; k<Nz; k++)
	for (j=0; j<Ny; j++)
	  for (i=0; i<Nx; i++)
	    W(h, k, j, i) = W(h, k-1, j, i)
	      + (
	       
		 ( U(h, k-1, j, i) - U(h, k-1, j, i+1) )
		 * (W[1].Value(h, k-1, j, i) - W[1].Value(h, k, j, i))
		 / (U[3].Value(h, k-1, j, i) - U[3].Value(h, k-1, j, i+1))

		 + ( V(h, k-1, j, i) - V(h, k-1, j+1, i) )
		 * (W[1].Value(h, k-1, j, i) - W[1].Value(h, k, j, i))
		 / (sin(pi*V[2].Value(h, k-1, j, i)/180.)
		    - sin(pi*V[2].Value(h, k-1, j+1, i)/180.))

		 ) / R_Earth;
  
  }

}  // namespace AtmoData.

#define ATMODATA_FILE_WINDS_CXX
#endif
