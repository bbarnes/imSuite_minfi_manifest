
#if HAVE_CONFIG_H
#include "config.h"
#endif

#include "imscore-run.h"

/* TO DO:
   - We should penalize the number of long runs of single residues.
   - We may also penalize long runs of two nucleotides. */

float
longest_run_score (const int length)
{
  if (length <= 4)
    return 1.00;
  else if (length <= 5)
    return 0.98;
  else if (length <= 6)
    return 0.95;
  else if (length <= 7)
    return 0.90;
  else if (length <= 8)
    return 0.85;
  else if (length <= 9)
    return 0.65;
  else if (length <= 10)
    return 0.55;
  else if (length <= 11)
    return 0.45;
  else if (length <= 12)
    return 0.35;
  else if (length <= 13)
    return 0.25;
  else
    return 0.10;
}

/*
  Local Variables:
  c-file-style: "gnu"
  End:
*/
