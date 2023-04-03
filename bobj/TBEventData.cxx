#include <ctime>
#include <sys/types.h>
#include <sys/stat.h>
#include "TBEventData.hxx"

ClassImp(TBEventData)

TBEventData::TBEventData() : TNamed("eventData", "eventData"){}

void TBEventData::update(time_t time) 
{
    evtime = time;
    struct tm *timeinfo;
    timeinfo = localtime(&evtime);
    year = timeinfo->tm_year;
    mon = timeinfo->tm_mon;
    day = timeinfo->tm_mday;
    hour= timeinfo->tm_hour;
    min = timeinfo->tm_min;
    sec = timeinfo->tm_sec;
    isdst = timeinfo->tm_isdst;
}