#include "anaRun.cc"
int main(int argc, char *argv[])
{
  cout << "executing " << argv[0] << endl;
  printf(" usage: ana  <run name>  <max files 0=all>  \n ");
  if(argc<2)
    exit(0);
  TString tag("run");
  Long64_t maxFiles = 0;
  if (argc > 1)
  {
    tag = TString(argv[1]);
  }
  if (argc > 2)
  {
    maxFiles = atoi(argv[2]);
  }

  printf(" starting anaRun %s %lld \n", tag.Data(), maxFiles);
  anaRun* r=new anaRun(tag);
  r->anaRunFile(tag,maxFiles);
  exit(0);
}
