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

}  // namespace AtmoData.


#define ATMODATA_FILE_FORMAT_HXX
#endif
