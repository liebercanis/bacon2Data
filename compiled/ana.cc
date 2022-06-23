#include "anaRun.cc"
int main(int argc, char* argv[])
{
  cout << "executing " << argv[0] << endl;
  if(argc<2) {
    printf(" usage: ana  <run name>   \n ");
    exit(0);
   }
  TString tag  = TString(argv[1]);

  printf(" starting anaRun %s \n", tag.Data());
  new anaRun(tag);
  exit(0);
}
