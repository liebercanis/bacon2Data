#include "TBFile.hxx"
ClassImp(TBFile)

TBFile::TBFile(TString name): TNamed(name,name)
{
  if (stat(name, &fileInfo) != 0)
  { // Use stat() to get the info
    std::cout << "Error: " << strerror(errno) << std::endl;
    return;
  }
  created=std::ctime(&fileInfo.st_ctime); // Creation time
  modified = std::ctime(&fileInfo.st_mtime);
  print();
}