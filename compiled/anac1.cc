#include "anaCRun.cc"
int main(int argc, char *argv[])
{
  cout << "executing " << argv[0] << endl;
  printf(" usage: ana  <run name>  <max entries  0=all>  \n ");
  if(argc<2)
    exit(0);
  TString tag("run");
  Long64_t maxEntries = 0;
  if (argc > 1)
  {
    tag = TString(argv[1]);
  }
  if (argc > 2)
  {
    maxEntries = atoi(argv[2]);
  }

  printf(" starting anaRun %s %lld \n", tag.Data(), maxEntries);
  anaCRun* r=new anaCRun(tag);
  r->anaCRunFile(tag,maxEntries);
  exit(0);
}
