
#include <assert.h>

#if HAVE_CONFIG_H
#include "config.h"
#endif

#include "imscore-cpg.h"

/* Ignore allele_type at this time */
float
underlying_cpg_score (const int cpg_count, const allele_type_t allele_type)
{
  assert (cpg_count >= 0);

  if (cpg_count == 0)
    return 1.00;
  else if (cpg_count == 1)
    return 0.99;
  else if (cpg_count == 2)
    return 0.98;
  else if (cpg_count == 3)
    return 0.97;
  else if (cpg_count == 4)
    return 0.50;
  else if (cpg_count == 5)
    return 0.25;
  else
    return 0.10;
}

/*
  Local Variables:
  c-file-style: "gnu"
  End:
*/
