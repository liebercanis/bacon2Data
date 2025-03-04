#include "TOyHit.hxx"
ClassImp(TOyHit)

    TOyHit::TOyHit() : TNamed("TOyHit", "TOyHit")
{
  clear();
}
void TOyHit::clear()
{
  event = 0;
  time = 0;
  val = 0;
};
