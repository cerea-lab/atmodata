#ifndef ATMODATA_FILE_ATMODATA_HXX

#include "SeldonData.hxx"
using namespace SeldonData;

template <class T>
void to_num(std::string s, T& num)
{
  std::istringstream str(s);
  str >> num;
}

template <class T>
T to_num(std::string s)
{
  T num;
  std::istringstream str(s);
  str >> num;
  return num;
}

#include "Kz.cxx"
#include "Winds.cxx"

#include "Format.cxx"

#define ATMODATA_FILE_ATMODATA_HXX
#endif
