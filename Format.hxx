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
