
#include <assert.h>
#include <string.h>
#include <stdio.h>

#if HAVE_CONFIG_H
#include "config.h"
#endif

#include "biolib.h"
#include "imdefs.h"
#include "imscore-self.h"

int
probe_self_any (const char *seq, const int seq_l)
{
  char r[MAX_PROBE_LENGTH], *l;
  int ret;
  dp_param_t a;
  dp_path_t o;

  if (seq == NULL || seq_l <= 0)
    return -1;

  assert (seq_l + 1 <= MAX_PROBE_LENGTH);
  nuseq_rc (seq, seq_l, r);

  init_dp_nt (&a);
  dp_global (seq, r, &a, &o);
  l = lcs (&o);

#ifdef DEBUG
  debug0 (__FILE__, __FUNCTION__, __LINE__, "Self complementarity");
  debug_s ("Sequence", seq);
  debug_s ("Reverse", r);
  dp_print_alignment_7 (&o, stderr);
  if (l)
    {
      debug_s ("LCS", l);
      debug_i ("LCS length", strlen (l));
    }
  else
    {
      debug_s ("LCS", l);
      debug_i ("LCS length", 0);
    }
#endif

  xfree (o.path[0]);
  xfree (o.path[1]);
  if (l == NULL)
    return 0;
  else
    {
      ret = strlen (l);
      xfree (l);
      return ret;
    }
}

/* Generally we do not want tax so-called self-complimentarity too much.

  This is a stem from a previously validated design.
  
  The decimal point is at the |

   0 | 0000000000000
   1 | 00000000000000000000000000000000000000000000000000000000000000000000+1169
   2 | 00000000000000000000000000000000000000000000000000000000000000000000+2872
   3 | 00000000000000000000000000000000000000000000000000000000000000000000+933
   4 | 00000000000000000000000000000000000000000000000000000000000000000000+560
   5 | 00000000000000000000000000000000000000000000000000000000000000000000
   6 | 00000000000000000000000000000000000000000000000000000000000000000000+54
   7 | 000
   8 | 000000000000000000000000000000000
   9 |
  10 | 000000
 */
float
probe_self_any_score (const float self_any)
{
  if (self_any <= 3.0)
    return 1.0;
  else if (self_any <= 4.0)
    return 0.99;
  else if (self_any <= 5.0)
    return 0.98;
  else if (self_any <= 6.0)
    return 0.96;
  else if (self_any <= 7.0)
    return 0.94;
  else if (self_any <= 8.0)
    return 0.90;
  else if (self_any <= 9.0)
    return 0.80;
  else
    return 0.10;
}

int
probe_self_end (const char *seq, const int seq_l)
{
  fprintf (stderr,
           "Self complimentarity END is neither used nor implemented.\n");
  return -1;
}

float
probe_self_end_score (const float self_end)
{
  fprintf (stderr,
           "Self complimentarity END is neither used nor implemented.\n");
  return 1.0;
}

/*
  Local Variables:
  c-file-style: "gnu"
  End:
*/
