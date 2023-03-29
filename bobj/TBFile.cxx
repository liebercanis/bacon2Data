#include <ctime>
#include <sys/types.h>
#include <sys/stat.h>
#include "TBFile.hxx"

ClassImp(TBFile)

TBFile::TBFile(): TNamed("tbfile","none")
{
}

TBFile::TBFile(TString name): TNamed(TString("tbfile"),name)
{
  update(name);
}

void TBFile::update(TString name){
  struct stat fileInfo;
  if (stat(name, &fileInfo) != 0)
  { // Use stat() to get the info
    std::cout << "Error: " << strerror(errno) << " using stat for file " << name << std::endl;
    return;
  }
  tsc = TTimeStamp(fileInfo.st_ctime);
  tsm = TTimeStamp(fileInfo.st_mtime);
  created = std::ctime(&fileInfo.st_ctime); // Creation time
  modified = std::ctime(&fileInfo.st_mtime);
  print();
}