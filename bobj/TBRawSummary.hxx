/**a
** MG, Oct 18, 2024
**/
#ifndef TBRAWSUMMARY_DEFINED
#define TBRAWSUMMARY_DEFINED
#include <iostream>
#include <string>
#include <TNamed.h>
using namespace std;
// class to store file info

class TBRawSummary : public TNamed
{
public:
  TBRawSummary();
  std::vector<double> mean;
  std::vector<double> sigma;

  void
  clear()
  {
    mean.clear();
    sigma.clear();
  }
  void print()
  {
    printf("TBRawSummary %lu channels \n", mean.size());
    if (mean.size() == sigma.size())
      for (unsigned long ic = 0; ic < mean.size(); ++ic)
        printf(" chan %lu mean %f sigma %f \n", ic, mean[ic], sigma[ic]);
  }
  ClassDef(TBRawSummary, 1)
};

#endif
