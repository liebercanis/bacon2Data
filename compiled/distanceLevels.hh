#ifndef DISTANCE_LEVEL_H
#define DISTANCE_LEVEL_H

enum
{
  NLEVELS = 5
};
bool geoVersionOld = true;

// distances in 3D radius of sipms
static double distanceLevel[NLEVELS];

static void setDistanceLevels(bool isOld = true)
{
  geoVersionOld = isOld;
  if (!geoVersionOld)
  {
    distanceLevel[0] = 1.062; // trigger
    distanceLevel[1] = 11.48;
    distanceLevel[2] = 21.41;
    distanceLevel[3] = 31.35;
    distanceLevel[4] = 42.70; // PMT
  }
  else
  {
    distanceLevel[0] = 1.251; // trigger
    distanceLevel[1] = 11.789;
    distanceLevel[2] = 21.723;
    distanceLevel[3] = 31.668;
    distanceLevel[4] = 42.017; // PMT
  }
}
#endif // DISTANCE_LEVEL_H

/*
distanceLevel[0] = 1.062; // 11.6;
distanceLevel[1] = 11.48; // 11.6;
distanceLevel[2] = 21.41; // 23.2;
distanceLevel[3] = 31.35; // 34.8;
distanceLevel[4] = 42.70; // 36.0;
*/

/* georgia doc https://www.overleaf.com/project/681e887700ae3302025fbaf6
  old
  r x y z
chan 9  1.251 0.795 -0.459 -0.851
chan 10 1.251 -0.795 -0.459 -0.851
chan 11 1.251 0.000 0.918 -0.851

chan 8 11.789 0.000 -2.722 -11.470
chan 7 11.789 2.357 1.361 -11.470
chan 6 11.789 -2.357 1.361 -11.470

chan 5 21.723 0.000 -3.566 -21.429
chan 4 21.723 3.089 1.783 -21.429
chan 3 21.723 -3.089 1.783 -21.429

chan 2 31.668 0.000 -3.657 -31.456
chan 1 31.668 3.167 1.828 -31.456
chan 0 31.668 -3.167 1.828 -31.456

PMT 43.017 0.000 0.000 -43.017

new
r x y z
chan 9 1.062 0.795 -0.459 -0.534
chan 10 1.062 -0.795 -0.459 -0.534
chan 11 1.062 0.000 0.918 -0.534

chan 8 11.480 0.000 -2.722 -11.153
chan 7 11.480 2.357 1.361 -11.153
chan 6 11.480 -2.357 1.361 -11.153

chan 5 21.411 0.000 -3.566 -21.112
chan 4 21.411 3.089 1.783 -21.112
chan 3 21.411 -3.089 1.783 -21.112

chan 2 31.353 0.000 -3.657 -31.139
chan 1 31.353 3.167 1.828 -31.139
chan 0 31.353 -3.167 1.828 -31.139

PMT 42.70 0.000 0.000 -42.700


*/