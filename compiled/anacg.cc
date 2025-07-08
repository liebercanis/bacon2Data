#include "anaCRunGamma.cc"

// static TBRun *theTBRun;
int main(int argc, char *argv[])
{
  std::cout << "executing gamma with anaCRunGamma derivative pulse finding " << argv[0] << std::endl;
  printf(" usage: ana  gamma <run name>  <max entries  0=all> <firstEvent default 0> \n ");
  if (argc < 2)
  {
    printf("... %i %s exit\n", argc, argv[0]);
    exit(1);
  }

  TString tag("run");
  Long64_t firstEntry = 0;
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

  // make output tree
  printf(" starting anaCRunGamma %s maxEntries %lld firstEntry %lld \n", tag.Data(), maxEntries, firstEntry);
  anaCRun *r = new anaCRun(tag);
  // r->setTBRun(theTBRun);
  int rc = r->anaCRunFile(tag, maxEntries, firstEntry);
  printf("... %s return code %i exit\n", argv[0], rc);
  if (rc == 0)
    exit(0);
  else
    exit(1);
}
