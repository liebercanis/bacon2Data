
int level(int ichan) {
   int ilevel = -1;
  if (ichan == 6 || ichan == 7 || ichan == 8)
    ilevel = 0;
  else if (ichan == 3 || ichan == 4 || ichan == 5)
    ilevel = 1;
  else if (ichan == 0 || ichan == 1 || ichan == 2)
    ilevel = 2;
  else if(ichan==12) 
    ilevel = 3;
  return ilevel;
}

double getDistance(int ilevel)
{
  if(ilevel<0) return 1.0;
  double distance[4];
  distance[0] = 11.6;
  distance[1] = 23.2;
  distance[2] = 34.8;;
  distance[3] = 36.0;
  return distance[ilevel];
}


double effLevel(int ilevel)
{
  /*
  **  Area of SiPMs is 6.0mm x 6.0mm

      Channels 6, 7, and 8 are at 11.6 cm
      from the source Channels 3, 4, and 5 are at 23.2 cm
      from the source Channels 0, 1, and 2 are at 34.8 cm from the source Channel 12 is at 36 cm from the source.
      */
  double a = pow(0.6, 2.);
  double b = 4.0 * TMath::Pi();
  double d2 = pow(getDistance(ilevel),2.);
  double e = a / b / d2;
  return e;
}

void effGeo(){


  for(int ichan=0; ichan<13; ++ ichan) {
    int ilevel = level(ichan);
    double d = getDistance(ilevel);
    double e = effLevel(ilevel);
    //printf("chan %i level %i dist %.3f e %.3E \n",ichan, ilevel, d, e);
    printf("effGeo[%i]=%f ; \n",ichan, e);
  }
}
