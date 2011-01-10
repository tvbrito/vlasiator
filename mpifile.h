#ifndef MPIWRITER_H
#define MPIWRITER_H

#include <mpi.h>

namespace MPIConv {
   // Template for converting a native C datatype into MPI datatype.
   template<typename T> inline MPI_Datatype MPItype();
   
   template<> inline MPI_Datatype MPItype<char>() {return MPI_CHAR;}
   template<> inline MPI_Datatype MPItype<float>() {return MPI_FLOAT;}
   template<> inline MPI_Datatype MPItype<double>() {return MPI_DOUBLE;}
   template<> inline MPI_Datatype MPItype<int>() {return MPI_INT;}
   template<> inline MPI_Datatype MPItype<long double>() {return MPI_LONG_DOUBLE;}
   template<> inline MPI_Datatype MPItype<long int>() {return MPI_LONG;}
   template<> inline MPI_Datatype MPItype<long long int>() {return MPI_LONG_LONG;}
   template<> inline MPI_Datatype MPItype<short int>() {return MPI_SHORT;}
   template<> inline MPI_Datatype MPItype<signed char>() {return MPI_SIGNED_CHAR;}
   template<> inline MPI_Datatype MPItype<unsigned char>() {return MPI_UNSIGNED_CHAR;}
   template<> inline MPI_Datatype MPItype<unsigned int>() {return MPI_UNSIGNED;}
   template<> inline MPI_Datatype MPItype<unsigned long int>() {return MPI_UNSIGNED_LONG;}
   template<> inline MPI_Datatype MPItype<unsigned long long int>() {return MPI_UNSIGNED_LONG_LONG;}
   template<> inline MPI_Datatype MPItype<unsigned short int>() {return MPI_UNSIGNED_SHORT;}
}

/**
 * Note that MPI must have been initialized prior to using this class.
 */
class MPIFile {
 public:
   MPIFile();
   ~MPIFile();
   
   bool open(MPI_Comm comm,const std::string& fname,MPI_Info file,const int& accessMode,const bool& deleteFile=false);
   bool close();
   MPI_Status getStatus() const {return MPIstatus;}
   bool resetPosition();
   
   template<typename T> MPIFile& operator<<(const T& value);
   template<typename T> MPIFile& operator>>(T& value);
   template<typename T> T getCount();
   template<typename T,typename INT> bool read(const INT& count,T* const array);
   template<typename T,typename INT> bool write(const INT& count,const T* const array);
   
 private:
   MPI_File fileptr;     /**< Pointer to file opened with MPI.*/
   MPI_Info MPIinfo;     /**< Variable for holding MPI file info.*/
   MPI_Status MPIstatus; /**< Variable for holding MPI status messages.*/
   void* ptr;            /**< Pointer used to write bytes to MPI files.*/
};

// ********************************************************
// ****** DEFINITIONS FOR MPIFILE TEMPLATE FUNCTIONS ******
// ********************************************************

/** Stream input operator for writing a primitive datatype into an MPI file. 
 * This function should work for all datatypes for which sizeof gives an accurate 
 * and portable size, i.e. do not use this function for writing structs or classes.
 * In such cases use MPIFile::write. It is mandatory to call MPIFile::resetPosition 
 * before attempting to access the file.
 * @param value The data to write.
 * @return A reference to MPIFile.
 * @see MPIFile::read.
 * @see MPIFile::resetPosition
 */
template<typename T> inline MPIFile& MPIFile::operator<<(const T& value) {
   ptr = static_cast<void*>(const_cast<T*>(&value));
   MPI_File_write_shared(fileptr,ptr,sizeof(T),MPI_BYTE,&MPIstatus);
   return *this;
}

/** Stream output operator for reading a primitive datatype from an MPI file.
 * This function should work for all datatypes for which sizeof gives an 
 * accurate and portable size, i.e. do not use this function for reading
 * structs or classes. Use MPIFIle::read instead. It is mandatory to call 
 * MPIFile::resetPosition before attampting to access the file.
 * @param value A variable in which the read value is written.
 * @return A reference to MPIFIle.
 * @see MPIFile::read
 * @see MPIFile::resetPosition
 */
template<typename T> inline MPIFile& MPIFile::operator>>(T& value) {
   ptr = static_cast<void*>(&value);
   MPI_File_read_shared(fileptr,ptr,sizeof(T),MPI_BYTE,&MPIstatus);
   return *this;
}

template<typename T> inline T MPIFile::getCount() {
   int count;
   MPI_Get_count(&MPIstatus,MPI_BYTE,&count);
   return static_cast<T>(count);
}

/** A function for reading an array of primitive datatypes from an MPI file. This function 
 * should be used to read non-trivial datatypes, such as structs and class data, 
 * by first reading a byte array containing the data and then typecasting it into correct 
 * variables. It is mandatory to call MPIFile::resetPosition before attempting to access 
 * the file.
 * @param count Number of elements to read.
 * @param array An array whose contents will be overwritten with the read data. This array 
 * must be allocated prior to calling this function.
 * @return If true, the data was read successfully.
 * @see MPIFile::resetPosition
 */
template<typename T,typename INT> inline bool MPIFile::read(const INT& count,T* const array) {
   if (MPI_File_read_shared(fileptr,array,sizeof(T)*count,MPI_BYTE,&MPIstatus) == MPI_SUCCESS) return true;
   return false;
}

/** A function for writing an array of primitive datatypes into an MPI file. This function
 * should be used to write non-trivial datatypes, such as structs and class data,
 * by first converting the data into a byte array. It is mandatory to call MPIFile::resetPosition 
 * before attempting to access the file.
 * @param count Number of elements in array.
 * @param array An array of primitive datatypes to write.
 * @return If true, the data was written successfully.
 * @see MPIFile::resetPosition
 */
template<typename T,typename INT> inline bool MPIFile::write(const INT& count,const T* const array) {
   //int mpiRank;
   //MPI_Comm_rank(MPI_COMM_WORLD,&mpiRank);
   //std::cerr << "\t Proc #" << mpiRank << " writing " << sizeof(T)*count << " bytes" << std::endl;
   
   ptr = static_cast<void*>(const_cast<T*>(array));
   //if (MPI_File_write_shared(fileptr,ptr,sizeof(T)*count,MPI_BYTE,&MPIstatus) == MPI_SUCCESS) return true;
   //return false;

   bool rvalue = true;
   if (MPI_File_write_shared(fileptr,ptr,sizeof(T)*count,MPI_BYTE,&MPIstatus) != MPI_SUCCESS) rvalue = false;
   /*
   MPI_Datatype cur_etype;
   MPI_Datatype cur_ftype;
   MPI_Offset cur_disp;
   char strarray[100];
   MPI_File_get_view(fileptr,&cur_disp,&cur_etype,&cur_ftype,strarray);
   std::cerr << "\t Current displacement after write is " << (unsigned long long int)cur_disp << std::endl;
   */
   //resetPosition();
   return rvalue;
}

#endif
