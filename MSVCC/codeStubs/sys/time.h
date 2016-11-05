
#ifndef TIME_STUB_H
#define TIME_STUB_H

#include <winsock.h>	//For timeval defn

struct timezone {
  int tz_minuteswest;
  int tz_dsttime;
};

int gettimeofday(struct timeval *tv, struct timezone *tz);


#endif