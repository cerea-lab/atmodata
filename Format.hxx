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


#ifndef ATMODATA_FILE_FORMAT_HXX

namespace AtmoData
{

  //! Input/ouput class to read files in CSV format (Airparif).
  class FormatCSV: public Format
  {

  protected:

  public:
    FormatCSV()  throw();
    ~FormatCSV()  throw();

    // Data.

    template<class TD, int N, class TG,
	     class TS, class TGS>
    void Read(string FileName, Data<TD, N, TG>& D,
	      Data<TS, 1, TGS>& S) const;
    template<class TD, int N, class TG,
	     class TS, class TGS>
    void Read(ifstream& FileStream, Data<TD, N, TG>& D,
	      Data<TS, 1, TGS>& S) const;

    // Array.

    template<class TA, int N,
	     class TS, class TGS>
    void Read(string FileName, Array<TA, N>& A,
	      Data<TS, 1, TGS>& S) const;
    template<class TA, int N,
	     class TS, class TGS>
    void Read(ifstream& FileStream, Array<TA, N>& A,
	      Data<TS, 1, TGS>& S) const;

  };  

  //! Input/ouput class to read files in binary format at ECMWF.
  template<class T>
  class FormatECMWF: public Format
  {

  protected:
    int date_;

  public:
    FormatECMWF()  throw();
    FormatECMWF(int date)  throw();
    ~FormatECMWF()  throw();

    void SetDate(int date);
    int GetDate() const;

    // Data.

    template<class TD, int N, class TG>
    void Read(string FileName, Data<TD, N, TG>& D) const;
    template<class TD, int N, class TG>
    void Read(ifstream& FileStream, Data<TD, N, TG>& D) const;

    // Array.

    template<class TA, int N>
    void Read(string FileName, Array<TA, N>& A) const;
    template<int N>
    void Read(ifstream& FileStream, Array<T, N>& A) const;
    template<class TA, int N>
    void Read(ifstream& FileStream, Array<TA, N>& A) const;

  };  

}  // namespace AtmoData.


#define ATMODATA_FILE_FORMAT_HXX
#endif
