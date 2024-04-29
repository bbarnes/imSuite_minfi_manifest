
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#if HAVE_CONFIG_H
#include "config.h"
#endif

#include "biolib.h"

#include "imdefs.h"
#include "imscore-tango.h"

/* #include "nmerscan.h" */

#define TANGO_WORDSIZE 11

static int *Freqs = NULL;
static int Wordsize = TANGO_WORDSIZE;

int
read_tango_11mer_data (const char *fn_11mer_b, const int nmer_wordsize)
{
  assert (nmer_wordsize == TANGO_WORDSIZE);     /* At this time */
  Wordsize = nmer_wordsize;

  if (fn_11mer_b == NULL)
    {
      fprintf (stderr, "Nmer-binary file name can not be NULL\n");
      exit (EXIT_FAILURE);
    }
  else
    {
      Freqs = nmer_read_fb (fn_11mer_b, nmer_wordsize);
      return 0;
    }
}

void
delete_tango_11mer_data ()
{
  xfree (Freqs);
}

int
tango_count (const char *seq, const int seq_l)
{
  int n_wins = MAX_PROBE_LENGTH, distr[MAX_PROBE_LENGTH];

  if (Freqs == NULL)
    {
      fprintf (stderr, "N-mer library is not read in.\n");
      return -1.0;
    }

  nmer_dist_dram (seq, seq_l, Freqs, Wordsize, &n_wins, &distr[0]);
  return sum_of_nmer_count (n_wins, distr);
}

float
tango_score (const int tango_count)
{
  float ret = -1.0;

  if (tango_count < 0)
    ret = -1.0;
  else if (tango_count == 0)
    ret = 1.0;
  else
    ret = 0.0;

  return ret;
}

/*
  Local Variables:
  c-file-style: "gnu"
  End:
*/
