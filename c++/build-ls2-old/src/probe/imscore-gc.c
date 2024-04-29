
#include <assert.h>
#include <math.h>

#if HAVE_CONFIG_H
#include "config.h"
#endif

#include "imdefs.h"
#include "imscore-gc.h"

float
gc_score (const float gc_percent, const allele_type_t allele_type)
{
  float ret;

  assert (allele_type == ALLELE_C || allele_type == ALLELE_T);
  assert (gc_percent >= 0.0);

  ret = 0.0;
  if (allele_type == ALLELE_C)  /* Methylated */
    {
      ret = exp (-0.00014 * pow (gc_percent - 39.00, 2.0)) + 0.001;
    }

  else if (allele_type == ALLELE_T)     /* UnMethylated */
    {
      ret = exp (-0.00010 * pow (gc_percent - 37.00, 2.0)) + 0.001;
    }

  if (ret > 1.0)
    ret = 1.0;
  else if (ret < 0.0)
    ret = 0.0;

  return ret;
}

/*
  Local Variables:
  c-file-style: "gnu"
  End:
*/
