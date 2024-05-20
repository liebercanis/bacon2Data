#include "anaCRun.cc"

static TBRun *theTBRun;
int main(int argc, char *argv[])
{
  cout << "executing " << argv[0] << endl;
  printf(" usage: ana  <run name>  <max entries  0=all> <firstEvent default 0> \n ");
  if (argc < 2)
    exit(0);
  TString tag("run");
  Long64_t firstEntry = 0;
  if (argc < 2)
    exit(0);
  Long64_t maxEntries = 0;
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
  theTBRun = new TBRun(tag);
  // make output tree
  printf(" starting anaRun %s maxEntries %lld firstEntry %lld \n", tag.Data(), maxEntries, firstEntry);
  anaCRun *r = new anaCRun(tag);
  //r->setTBRun(theTBRun);
  r->anaCRunFile(tag, maxEntries, firstEntry);
  exit(0);
}
