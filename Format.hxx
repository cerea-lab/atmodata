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

  // For MM5 sub-headers.
  class MM5SubHeader;

  //! Input/ouput class to read MM5 files (version 3).
  class FormatMM5: public Format
  {

  protected:

  public:
    FormatMM5()  throw();
    ~FormatMM5()  throw();

    // Flag.
    int ReadFlag(ifstream& FileStream) const;

    // Big header.
    void ReadBigHeader(string FileName,
		       Array<int, 2>& BHI, Array<float, 2>& BHR,
		       Array<string, 2>& BHIC, Array<string, 2>& BHRC) const;
    void ReadBigHeader(ifstream& FileStream,
		       Array<int, 2>& BHI, Array<float, 2>& BHR,
		       Array<string, 2>& BHIC, Array<string, 2>& BHRC) const;
    void ReadBigHeader(ifstream& FileStream) const;

    // Sub-header.
    void ReadSubHeader(ifstream& FileStream, MM5SubHeader& SH) const;
    void ReadSubHeader(ifstream& FileStream) const;

    // Field.
    template <int N, class TG>
    void ReadWholeField(string FileName, string FieldName, 
			Data<float, N, TG>& A) const;
    template <int N>
    void ReadWholeField(string FileName, string FieldName, 
			Array<float, N>& A) const;
    template <int N, class TG>
    void ReadWholeField(ifstream& FileStream, string FieldName, 
			Data<float, N, TG>& A) const;
    template <int N>
    void ReadWholeField(ifstream& FileStream, string FieldName, 
			Array<float, N>& A) const;
    template <int N>
    void ReadField(ifstream& FileStream, bool cross, Array<float, N>& A) const;
    template <int N, class TG>
    void ReadField(ifstream& FileStream, Data<float, N, TG>& A) const;
    template <int N>
    void ReadField(ifstream& FileStream, MM5SubHeader& SH, Array<float, N>& A) const;
    template <int N>
    void ReadField(ifstream& FileStream, Array<float, N>& A) const;
    void ReadField(ifstream& FileStream) const;

  };

  //! For MM5 sub-headers.
  class MM5SubHeader
  {

  public:
    int ndim;
    Array<int, 1> start_index;
    Array<int, 1> end_index;
    float xtime;
    string staggering;
    string ordering;
    string current_date;
    string name;
    string unit;
    string description;

  public:
    MM5SubHeader()  throw();
    MM5SubHeader(const MM5SubHeader&)  throw();
    ~MM5SubHeader()  throw();

    void Init();
    MM5SubHeader& operator=(MM5SubHeader&);

  };

}  // namespace AtmoData.


#define ATMODATA_FILE_FORMAT_HXX
#endif
