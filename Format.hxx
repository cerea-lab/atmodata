#ifndef ATMODATA_FILE_FORMAT_HXX

namespace AtmoData
{

  //! Input/ouput class to read files in CSV format.
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

}  // namespace AtmoData.


#define ATMODATA_FILE_FORMAT_HXX
#endif
