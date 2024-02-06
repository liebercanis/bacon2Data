#include "anaRun.cc"
int main(int argc, char *argv[])
{
  cout << "executing " << argv[0] << endl;
  printf(" usage: ana  <run name>  <max entries  0=all> <firstEvent default 0> \n ");
  if(argc<2)
    exit(0);
  TString tag("run");
  Long64_t maxEntries = 0;
  Long64_t firstEntry = 0;
  if (argc > 1)
  {
    tag = TString(argv[1]);
  }
  if (argc > 2)
  {
    maxEntries = atoi(argv[2]);
  }
  if (argc > 3)
  {
    firstEntry = atoi(argv[3]);
  }

  printf(" starting anaRun %s %lld \n", tag.Data(), maxEntries, firstEntry);
  anaRun* r=new anaRun(tag);
  r->anaRunFile(tag,maxEntries);
  exit(0);
}
