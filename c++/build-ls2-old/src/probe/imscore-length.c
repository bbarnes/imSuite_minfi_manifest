
#if HAVE_CONFIG_H
#include "config.h"
#endif

#include "imscore-length.h"

float
probe_length_score (const int probe_length)
{
  if (probe_length == 50)
    return 1.0;
  else
    return 0.0;
}

/*
  Local Variables:
  c-file-style: "gnu"
  End:
*/
