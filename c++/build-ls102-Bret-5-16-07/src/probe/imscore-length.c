
#if HAVE_CONFIG_H
#include "config.h"
#endif

#include "imdefs.h"
#include "imscore-length.h"

/* In Infinium Methylation design, the probe_length is fixed.  The
   zero score for the nonconforming probes is a guard ensuing valid
   designs */
float
probe_length_score (const int probe_length)
{
  if (probe_length == INF_PROBE_LENGTH)
    return 1.0;
  else
    return 0.0;
}

/*
  Local Variables:
  c-file-style: "gnu"
  End:
*/
