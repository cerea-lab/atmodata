#ifndef ATMODATA_FILE_WINDS_HXX

namespace AtmoData
{

  template<class TU, class TV, class TW, class TG>
  void GetVerticalWind(Data<TU, 4, TG>& U, Data<TV, 4, TG>& V, Data<TW, 4, TG>& W);

}  // namespace AtmoData.

#define ATMODATA_FILE_WINDS_HXX
#endif
