#include <cstdlib>
#include <iostream>

#include "mpifile.h"

using namespace std;

/** Open a file for parallel writing with MPI.
 * @param comm The communicator describing the processes which will access the file. Use 
 * MPI_COMM_WORLD if every process needs to accesse the file. 
 * @param fname The name of the file.
 * @param accessMode Description of the mode in which the file is opened. See MPI specification.
 * @return If true, the file was opened successfully.
 */
bool MPIFile::open(MPI_Comm comm,const std::string& fname,MPI_Info info,const int& accessMode,const bool& deleteFile) {
   MPIinfo = info;
   if (deleteFile == true) MPI_File_delete(const_cast<char*>(fname.c_str()),MPI_INFO_NULL);
   if (MPI_File_open(comm,const_cast<char*>(fname.c_str()),accessMode,MPIinfo,&fileptr) != MPI_SUCCESS) {
      fileptr = MPI_FILE_NULL;
      return false;
   }
   return true;
}

/** Close a file which has been previously opened with MPIFile::open.
 * @return If true, the file was closed without errors.
 * @see MPIFile::open
 */
bool MPIFile::close() {
   if (MPI_File_close(&fileptr) == MPI_SUCCESS) return true;
   return false;
}

/** Reset the file position pointer and MPI file view. This function must 
 * be called before writing to the file. Please see MPI specification for 
 * more information about file view.
 * @return If true, the file view was set without errors.
 */
bool MPIFile::resetPosition() {
   if (MPI_File_set_view(fileptr,MPI_DISPLACEMENT_CURRENT,MPI_BYTE,MPI_BYTE,const_cast<char*>("native"),MPIinfo) != MPI_SUCCESS) {
      std::cerr << "MPIwriter ERROR: Failed to set file view for process 0" << std::endl;
      return false;
   }
   return true;
}

// ****************************************
// ****** CONSTRUCTORS & DESTRUCTORS ******
// ****************************************

MPIFile::MPIFile() {
   
}

MPIFile::~MPIFile() {
   
}




