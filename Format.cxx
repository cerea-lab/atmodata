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

}  // namespace AtmoData.


#define ATMODATA_FILE_FORMAT_HXX
#endif
