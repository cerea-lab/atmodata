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
