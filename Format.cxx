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


#ifndef ATMODATA_FILE_FORMAT_HXX

#include "Format.hxx"

namespace AtmoData
{

  ///////////////
  // FORMATCVS //
  ///////////////

  //! Default constructor.
  FormatCSV::FormatCSV()  throw()
  {
  }

  //! Destructor.
  FormatCSV::~FormatCSV()  throw()
  {
  }

  /********/
  /* Data */
  /********/
 
  //! Reads a file in "CSV" format (Airparif).
  template<class TD, int N, class TG, class TS, class TGS>
  void FormatCSV::Read(string FileName, Data<TD, N, TG>& D,
		       Data<TS, 1, TGS>& S) const
  {

    this->Read(FileName, D.GetArray(), S);

  }

  //! Reads a file in "CSV" format (Airparif).
  template<class TD, int N, class TG, class TS, class TGS>
  void FormatCSV::Read(ifstream& FileStream, Data<TD, N, TG>& D,
		       Data<TS, 1, TGS>& S) const
  {

    this->Read(FileStream, D.GetArray(), S);

  }

  /*********/
  /* Array */
  /*********/

  //! Reads a file in "CSV" format (Airparif).
  template<class TA, int N, class TS, class TGS>
  void FormatCSV::Read(string FileName, Array<TA, N>& A,
		       Data<TS, 1, TGS>& S) const
  {

    ifstream FileStream;
    FileStream.open(FileName.c_str(), ifstream::in);
    FileStream.flags(fstream::skipws);

#ifdef DEBUG_SELDONDATA_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("FormatCSV::Read(string FileName, Array<TA, N>& A)",
		    "Unable to open file \"" + FileName + "\".");
#endif

    this->Read(FileStream, A, S);

    FileStream.close();

  }

  //! Reads a file in "CSV" format (Airparif).
  template<class TA, int N, class TS, class TGS>
  void FormatCSV::Read(ifstream& FileStream, Array<TA, N>& A,
		       Data<TS, 1, TGS>& S) const
  {

#ifdef DEBUG_SELDONDATA_IO
    // Checks if the file ready.
    if (!FileStream.good())
      throw IOError("FormatCSV::Read(ifstream& FileStream, Array<TA, N>& A)",
		    "File is not ready.");
#endif

    int i, j, k;

    int nb_elements = A.numElements();
    int nb_stations = A.extent(0);
    int nb_meas = nb_elements / nb_stations;

    string temp_line, temp, temp_str;
    
    getline(FileStream, temp_line);
    getline(FileStream, temp_line);

    Array<int, 1> Index(10), Length(10);
    for (j=0; j<10; j++)
      {
	Index(j) = 0;
	Length(j) = A.extent(j);
      }

    i = 0;
    j = N-1;
    while ( (FileStream.good()) && (i<nb_elements) )
      {

	// Station name.
	getline(FileStream, temp_str, ',');

	k = 0;
	while ( (k<nb_stations) && (strstr(S(k).c_str(), temp_str.c_str())) )
	  k++;

	if (k==nb_stations)
	  {
	    cout << temp_str + ": station not found!" << endl;
	    throw "Error";
	  }

	// Pollutant.
	getline(FileStream, temp_str, ',');
	// Unit.
	getline(FileStream, temp_str, ',');
	// Date.
	getline(FileStream, temp_str, ',');

	k = 0;
	while ( (FileStream.good()) && (k<nb_meas) )
	  {
	    // Reads useless fields.
	    if ( (k%24==0) && (k!=0) )
	      {
		// Average.
		getline(FileStream, temp_str);
		// Station name.
		getline(FileStream, temp_str, ',');
		// Pollutant.
		getline(FileStream, temp_str, ',');
		// Unit.
		getline(FileStream, temp_str, ',');
		// Date.
		getline(FileStream, temp_str, ',');
	      }
	    getline(FileStream, temp_str, ',');

	    A(Index(0), Index(1), Index(2),
	      Index(3), Index(4), Index(5),
	      Index(6), Index(7), Index(8),
	      Index(9)) = to_num<TA>(temp_str);

	    j = N-1;
	    while ( (j>=0) && (Index(j)==Length(j)-1) )
	      {
		Index(j) = 0;
		j--;
	      }
	    
	    if (j!=-1)
	      Index(j)++;
	    
	    k++;
	  }

      }

  }


  /////////////////
  // FORMATECMWF //
  /////////////////

  //! Default constructor.
  template<class T>
  FormatECMWF<T>::FormatECMWF()  throw()
  {
    date_ = -1;
  }

  //! Constructor.
  template<class T>
  FormatECMWF<T>::FormatECMWF(int date)  throw()
  {
    date_ = date;
  }

  //! Destructor.
  template<class T>
  FormatECMWF<T>::~FormatECMWF()  throw()
  {
  }

  //! Sets the date.
  /*!
    \param date date.
  */
  template<class T>
  void FormatECMWF<T>::SetDate(int date)
  {
    date_ = date;
  }

  //! Get the date.
  /*!
    \return The date.
  */
  template<class T>
  int FormatECMWF<T>::GetDate() const
  {
    return date_;
  }

  /********/
  /* Data */
  /********/
 
  //! Reads a file in "ECMWF" format.
  template<class T>
  template<class TD, int N, class TG>
  void FormatECMWF<T>::Read(string FileName, Data<TD, N, TG>& D) const
  {

    this->Read(FileName, D.GetArray());

  }

  //! Reads a file in "ECMWF" format.
  template<class T>
  template<class TD, int N, class TG>
  void FormatECMWF<T>::Read(ifstream& FileStream, Data<TD, N, TG>& D) const
  {

    this->Read(FileStream, D.GetArray());

  }

  /*********/
  /* Array */
  /*********/

  //! Reads a file in "ECMWF" format.
  template<class T>
  template<class TA, int N>
  void FormatECMWF<T>::Read(string FileName, Array<TA, N>& A) const
  {

    ifstream FileStream;
    FileStream.open(FileName.c_str(), ifstream::in);
    FileStream.flags(fstream::skipws);

#ifdef DEBUG_SELDONDATA_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("FormatECMWF::Read(string FileName, Array<T, N>& A)",
		    "Unable to open file \"" + FileName + "\".");
#endif

    this->Read(FileStream, A);

    FileStream.close();

  }

  //! Reads a file in "ECMWF" format.
  template<class T>
  template<int N>
  void FormatECMWF<T>::Read(ifstream& FileStream, Array<T, N>& A) const
  {

    int nb_elements = A.numElements();
    unsigned long data_size = nb_elements * sizeof(T);
    int Nt = A.extent(0);
    int rec_length = data_size / Nt;

    T* data = A.data();

#ifdef DEBUG_SELDONDATA_IO

    // Checks if the file ready.
    if (!FileStream.good())
      throw IOError("FormatECMWF::Read(ifstream& FileStream, Array<T, N>& A)",
		    "File is not ready.");

    // Checks records length.
    streampos position;
    position = FileStream.tellg();

    int file_rec_length;
    FileStream.read(reinterpret_cast<char*>(&file_rec_length), 4);
    if (rec_length!=(file_rec_length-4))
      throw IOError("FormatECMWF<T>::Read(ifstream& FileStream, Array<T, N>& A)",
		    "Record length (as in file) is " + to_str(file_rec_length)
		    + " byte(s)," + " but data record length is "
		    + to_str(rec_length) + " byte(s) long.");

    FileStream.seekg(position);

#endif

    int i = 0;

    int date;
    bool reading = (date_==-1);

    while ( (!reading) && (FileStream.good()) )
      {
	// Record length.
	FileStream.read(reinterpret_cast<char*>(&date), 4);
	// Date.
	FileStream.read(reinterpret_cast<char*>(&date), 4);

	// Data.
	FileStream.read(reinterpret_cast<char*>(data), rec_length);
	reading = (date==date_);

	// Record length.
	FileStream.read(reinterpret_cast<char*>(&date), 4);

	i = 1;
      }
    
#ifdef DEBUG_SELDONDATA_IO

    // Checks if all was read.
    if (!reading)
      throw IOError("FormatECMWF::Read(ifstream& FileStream, Array<T, N>& A)",
		    "The date was not found.");

    // Checks file length.
    position = FileStream.tellg();
    FileStream.seekg(0, ios::end);
    unsigned long file_size = FileStream.tellg() - position;

    if (data_size>(file_size-12*(Nt-i)+i*rec_length))
      {
	throw IOError("FormatBinary<T>::Read(ifstream& FileStream, Array<T, N>& A)",
		      "Unable to read " + to_str(data_size) + " byte(s)." +
		      " The input stream is only "
		      + to_str((file_size/(rec_length+12)+i)*rec_length) + " byte(s) long.");
      }
    FileStream.seekg(position);

#endif

    while ( (i<Nt) && (FileStream.good()) )
      {
	// Record length.
	FileStream.read(reinterpret_cast<char*>(&date), 4);
	// Date.
	FileStream.read(reinterpret_cast<char*>(&date), 4);

	// Data.
	FileStream.read(reinterpret_cast<char*>(data) + rec_length * i,
			rec_length);

	// Record length.
	FileStream.read(reinterpret_cast<char*>(&date), 4);
	reading = (date==date_);

	i++;
      }
    
  }

  //! Reads a file in "ECMWF" format.
  template<class T>
  template<class TA, int N>
  void FormatECMWF<T>::Read(ifstream& FileStream, Array<TA, N>& A) const
  {

    int nb_elements = A.numElements();
    unsigned long data_size = nb_elements * sizeof(TA);
    int Nt = A.extent(0);
    int rec_length = data_size / Nt;
    int length = rec_length/sizeof(TA);

    TA* data_output = A.data();

    int input_rec_length = length*sizeof(T);
    T* data = new T[length];

#ifdef DEBUG_SELDONDATA_IO

    // Checks if the file ready.
    if (!FileStream.good())
      throw IOError("FormatECMWF::Read(ifstream& FileStream, Array<TA, N>& A)",
		    "File is not ready.");

    // Checks records length.
    streampos position;
    position = FileStream.tellg();

    int file_rec_length;
    FileStream.read(reinterpret_cast<char*>(&file_rec_length), 4);
    if (input_rec_length!=(file_rec_length-4))
      throw IOError("FormatECMWF<T>::Read(ifstream& FileStream, Array<TA, N>& A)",
		    "Record length (as in file) is "+ to_str(file_rec_length/sizeof(T))
		    + " element(s)," + " but data record length is "
		    + to_str(length) + " element(s) long.");

    FileStream.seekg(position);

#endif

    int i = 0;
    int j;

    int date;
    bool reading = (date_==-1);

    while ( (!reading) && (FileStream.good()) )
      {
	// Record length.
	FileStream.read(reinterpret_cast<char*>(&date), 4);
	// Date.
	FileStream.read(reinterpret_cast<char*>(&date), 4);

	// Data.
	FileStream.read(reinterpret_cast<char*>(data), input_rec_length);
	reading = (date==date_);

	if (reading)
	  for (j=0; j<length; j++)
	    data_output[j] = data[j];

	// Record length.
	FileStream.read(reinterpret_cast<char*>(&date), 4);

	i = 1;
      }
    
#ifdef DEBUG_SELDONDATA_IO

    // Checks if all was read.
    if (!reading)
      throw IOError("FormatECMWF::Read(ifstream& FileStream, Array<TA, N>& A)",
		    "The date was not found.");

    // Checks file length.
    position = FileStream.tellg();
    FileStream.seekg(0, ios::end);
    unsigned long file_size = FileStream.tellg() - position;

    if (data_size>((file_size-12*(Nt-i))/sizeof(T)*sizeof(TA)+i*rec_length))
      {
	throw IOError("FormatBinary<T>::Read(ifstream& FileStream, Array<T, N>& A)",
		      "Unable to read " + to_str(nb_elements) + " elements(s)."
		      + " The input stream is only "
		      + to_str((file_size/(12+input_rec_length)+i)*length) + " elements(s) long.");
      }
    FileStream.seekg(position);

#endif

    while ( (i<Nt) && (FileStream.good()) )
      {
	// Record length.
	FileStream.read(reinterpret_cast<char*>(&date), 4);
	// Date.
	FileStream.read(reinterpret_cast<char*>(&date), 4);

	// Data.
	FileStream.read(reinterpret_cast<char*>(data), input_rec_length);
	for (j=0; j<length; j++)
	  data_output[j + i*length] = data[j];

	// Record length.
	FileStream.read(reinterpret_cast<char*>(&date), 4);
	reading = (date==date_);

	i++;
      }
    
  }


  ///////////////
  // FORMATMM5 //
  ///////////////

  //! Default constructor.
  FormatMM5::FormatMM5()  throw()
  {
  }

  //! Destructor.
  FormatMM5::~FormatMM5()  throw()
  {
  }

  /********/
  /* Flag */
  /********/

  //! Reads a flag.
  int FormatMM5::ReadFlag(ifstream& FileStream) const
  {

#ifdef DEBUG_SELDONDATA_IO
    // Checks if the file ready.
    if (!FileStream.good())
      throw IOError("FormatMM5<T>::ReadFlag(ifstream& FileStream)",
		    "File is not ready.");
#endif

    int length, flag;

    FileStream.read(reinterpret_cast<char*>(&length), 4);
    FileStream.read(reinterpret_cast<char*>(&flag), 4);
    FileStream.read(reinterpret_cast<char*>(&length), 4);

#ifdef DEBUG_SELDONDATA_IO
    // Checks if the flag was read.
    if (!FileStream.good())
      throw IOError("FormatMM5<T>::ReadFlag(ifstream& FileStream)",
		    "Unable to read the flag.");
#endif

    return swap(flag);

  }

  /*************/
  /* BigHeader */
  /*************/

  //! Reads big header.
  void FormatMM5::ReadBigHeader(string FileName,
				Array<int, 2>& BHI, Array<float, 2>& BHR,
				Array<string, 2>& BHIC, Array<string, 2>& BHRC) const
  {

    ifstream FileStream;
    FileStream.open(FileName.c_str(), ifstream::binary);

#ifdef DEBUG_SELDONDATA_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("FormatMM5<T>::ReadBigHeader(string FileName, Array<int, 2>&, Array<float, 2>&, Array<string, 2>&, Array<string, 2>&)",
		    "Unable to open file \"" + FileName + "\".");
#endif

    this->ReadFlag(FileStream);
    this->ReadBigHeader(FileStream, BHI, BHR, BHIC, BHRC);

    FileStream.close();

  }

  //! Reads big header.
  void FormatMM5::ReadBigHeader(ifstream& FileStream,
				Array<int, 2>& BHI, Array<float, 2>& BHR,
				Array<string, 2>& BHIC, Array<string, 2>& BHRC) const
  {

#ifdef DEBUG_SELDONDATA_IO
    // Checks if the file is ready.
    if (!FileStream.good())
      throw IOError("FormatMM5<T>::ReadBigHeader(ifstream& FileStream, Array<int, 2>&, Array<float, 2>&, Array<string, 2>&, Array<string, 2>&)",
		    "File is not ready.");
#endif

    int i, j;
    int length;
    FileStream.read(reinterpret_cast<char*>(&length), 4);

    BHI.resize(20, 50);
    BHR.resize(20, 20);
    BHIC.resize(50, 20);
    BHRC.resize(20, 20);

    FileStream.read(reinterpret_cast<char*>(BHI.data()), 1000 * sizeof(int));
    for (i=0; i<20; i++)
      for (j=0; j<50; j++)
	swap(BHI(i, j));

    FileStream.read(reinterpret_cast<char*>(BHR.data()), 400 * sizeof(float));
    for (i=0; i<20; i++)
      for (j=0; j<20; j++)
	swap(BHR(i, j));

    for (i=0; i<20; i++)
      for (j=0; j<50; j++)
	{
	  BHIC(i, j).resize(80);
	  FileStream.read(const_cast<char*>(BHIC(i, j).c_str()), 80 * sizeof(char));
	}

    for (i=0; i<20; i++)
      for (j=0; j<20; j++)
	{
	  BHRC(i, j).resize(80);
	  FileStream.read(const_cast<char*>(BHRC(i, j).c_str()), 80 * sizeof(char));
	}

    FileStream.read(reinterpret_cast<char*>(&length), 4);

#ifdef DEBUG_SELDONDATA_IO
    // Checks if all was read.
    if (!FileStream.good())
      throw IOError("FormatMM5<T>::ReadBigHeader(ifstream& FileStream, Array<int, 2>&, Array<float, 2>&, Array<string, 2>&, Array<string, 2>&)",
		    "Unable to read the big header.");
#endif

  }

  //! Reads big header.
  void FormatMM5::ReadBigHeader(ifstream& FileStream) const
  {

    Array<int, 2> BHI;
    Array<float, 2> BHR;
    Array<string, 2> BHIC;
    Array<string, 2> BHRC;
    
    this->ReadBigHeader(FileStream, BHI, BHR, BHIC, BHRC);

  }

  //! Reads sub-header.
  void FormatMM5::ReadSubHeader(ifstream& FileStream, MM5SubHeader& SH) const
  {

#ifdef DEBUG_SELDONDATA_IO
    // Checks if the file is ready.
    if (!FileStream.good())
      throw IOError("FormatMM5<T>::ReadSubHeader(ifstream& FileStream, MM5SubHeader& SH)",
		    "File is not ready.");
#endif

    int length;
    FileStream.read(reinterpret_cast<char*>(&length), 4);

    SH.Init();

    FileStream.read(reinterpret_cast<char*>(&SH.ndim), sizeof(int));
    swap(SH.ndim);
    FileStream.read(reinterpret_cast<char*>(SH.start_index.data()), 4 * sizeof(int));
    swap(SH.start_index);
    FileStream.read(reinterpret_cast<char*>(SH.end_index.data()), 4 * sizeof(int));
    swap(SH.end_index);
    FileStream.read(reinterpret_cast<char*>(&SH.xtime), sizeof(float));
    swap(SH.xtime);
    FileStream.read(const_cast<char*>(SH.staggering.c_str()), 4 * sizeof(char));
    FileStream.read(const_cast<char*>(SH.ordering.c_str()), 4 * sizeof(char));
    FileStream.read(const_cast<char*>(SH.current_date.c_str()), 24 * sizeof(char));
    FileStream.read(const_cast<char*>(SH.name.c_str()), 9 * sizeof(char));
    FileStream.read(const_cast<char*>(SH.unit.c_str()), 25 * sizeof(char));
    FileStream.read(const_cast<char*>(SH.description.c_str()), 46 * sizeof(char));

    FileStream.read(reinterpret_cast<char*>(&length), 4);

#ifdef DEBUG_SELDONDATA_IO
    // Checks if all was read.
    if (!FileStream.good())
      throw IOError("FormatMM5<T>::ReadSubHeader(ifstream& FileStream, MM5SubHeader& SH)",
		    "Unable to read the sub-header.");
#endif

  }

  //! Reads sub-header.
  void FormatMM5::ReadSubHeader(ifstream& FileStream) const
  {

    MM5SubHeader SH;
    this->ReadSubHeader(FileStream, SH);

  }

  /*********/
  /* Field */
  /*********/

  //! Reads the field.
  template <int N, class TG>
  void FormatMM5::ReadField(ifstream& FileStream, Data<float, N, TG>& D) const
  {
    this->ReadField(FileStream, D.GetArray());
  }

  //! Reads the field.
  template <int N>
  void FormatMM5::ReadField(ifstream& FileStream, Array<float, N>& A) const
  {

    unsigned long data_size = A.numElements() * sizeof(float);
    float* data = A.data();

#ifdef DEBUG_SELDONDATA_IO
    // Checks if the file is ready.
    if (!FileStream.good())
      throw IOError("FormatMM5<T>::ReadField(ifstream& FileStream, Array<float, N>&)",
		    "File is not ready.");

    // Checks file length.
    streampos position;
    position = FileStream.tellg();
    FileStream.seekg(0, ios::end);
    unsigned long file_size = FileStream.tellg() - position;

    if (data_size + 8 > file_size)
      throw IOError("FormatMM5<T>::ReadField(ifstream& FileStream, Array<float, N>& A)",
		    "Unable to read (" + to_str(data_size) + " + 8) byte(s)." +
		    " The input stream is only " + to_str(file_size) + " byte(s) long.");

    FileStream.seekg(position);
#endif

    int length;
    FileStream.read(reinterpret_cast<char*>(&length), 4);

    FileStream.read(reinterpret_cast<char*>(data), data_size);
    swap(A);

    FileStream.read(reinterpret_cast<char*>(&length), 4);

#ifdef DEBUG_SELDONDATA_IO
    // Checks if all was read.
    if (!FileStream.good())
      throw IOError("FormatMM5<T>::ReadField(ifstream& FileStream, Array<float, N>&)",
		    "Unable to read the field.");
#endif

  }

  //! Reads the field.
  void FormatMM5::ReadField(ifstream& FileStream) const
  {

#ifdef DEBUG_SELDONDATA_IO
    // Checks if the file is ready.
    if (!FileStream.good())
      throw IOError("FormatMM5<T>::ReadField(ifstream& FileStream)",
		    "File is not ready.");
#endif

    unsigned long data_size;
    FileStream.read(reinterpret_cast<char*>(&data_size), 4);
    swap(data_size);

#ifdef DEBUG_SELDONDATA_IO
    // Checks if the file is ready.
    if (!FileStream.good())
      throw IOError("FormatMM5<T>::ReadField(ifstream& FileStream)",
		    "File is not ready.");

    // Checks file length.
    streampos position;
    position = FileStream.tellg();
    FileStream.seekg(0, ios::end);
    unsigned long file_size = FileStream.tellg() - position;

    if (data_size + 4 > file_size)
      throw IOError("FormatMM5<T>::ReadField(ifstream& FileStream)",
		    "Unable to read (" + to_str(data_size) + " + 4) byte(s)." +
		    " The input stream is only " + to_str(file_size) + " byte(s) long.");

    FileStream.seekg(position);
#endif

    FileStream.seekg(data_size, ios::cur);

    unsigned long length;
    FileStream.read(reinterpret_cast<char*>(&length), 4);
    swap(length);

#ifdef DEBUG_SELDONDATA_IO
    // Checks if all was read.
    if (!FileStream.good())
      throw IOError("FormatMM5<T>::ReadField(ifstream& FileStream)",
		    "Unable to read the field.");
    // Checks if record length is the same as 'data_size'.
    if (length!=data_size)
      throw IOError("FormatMM5<T>::ReadField(ifstream& FileStream)",
		    "Unable to get the field size.");
#endif

  }


  //////////////////
  // MM5SUBHEADER //
  //////////////////

  //! Default constructor.
  MM5SubHeader::MM5SubHeader()  throw()
  {
  }

  //! Copy constructor.
  MM5SubHeader::MM5SubHeader(const MM5SubHeader& SH)  throw():
    ndim(SH.ndim),
    start_index(SH.start_index.copy()),
    end_index(SH.end_index.copy()),
    xtime(SH.xtime),
    staggering(SH.staggering),
    ordering(SH.ordering),
    current_date(SH.current_date),
    name(SH.name),
    unit(SH.unit),
    description(SH.description)
  {
  }

  //! Destructor.
  MM5SubHeader::~MM5SubHeader()  throw()
  {
  }

  //! Initializes sizes.
  void MM5SubHeader::Init()
  {
    start_index.resize(4);
    end_index.resize(4);
    staggering.resize(4);
    ordering.resize(4);
    current_date.resize(24);
    name.resize(9);
    unit.resize(25);
    description.resize(46);
  }

  //! Copies a sub-header.
  MM5SubHeader& MM5SubHeader::operator=(MM5SubHeader& SH)
  {
    this->ndim = SH.ndim;
    this->start_index = SH.start_index.copy();
    this->end_index = SH.end_index.copy();
    this->xtime = SH.xtime;
    this->staggering = SH.staggering;
    this->ordering = SH.ordering;
    this->current_date = SH.current_date;
    this->name = SH.name;
    this->unit = SH.unit;
    this->description = SH.description;
  }


}  // namespace AtmoData.


#define ATMODATA_FILE_FORMAT_HXX
#endif
