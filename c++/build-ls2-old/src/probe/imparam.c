
/* by Lixin Zhou */

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <unistd.h>
#include <sys/utsname.h>

#if HAVE_CONFIG_H
#include "config.h"
#endif

#include "biolib.h"

#include "imdefs.h"
#include "imparam.h"

/* These may (and perhaps should) go to config.h or configure.ac */
#define LEN_MIN_DEFAULT 50
#define LEN_MAX_DEFAULT 50
#define GC_MIN_DEFAULT 20
#define GC_MAX_DEFAULT 80
#define MIN_SCORE_TO_PRINT_DEFAULT 0.0

/* This structure is meant to be reusable for all the sequences in a
   single batch of probe design.  Some fields are to be updated for
   each sequence. */
probe_param_t *
probe_param_new ()
{
  probe_param_t *ppt;

  ppt = (probe_param_t *) xmalloc (sizeof (probe_param_t));

  ppt->len_min = LEN_MIN_DEFAULT;
  ppt->len_max = LEN_MAX_DEFAULT;
  ppt->gc_min = GC_MIN_DEFAULT;
  ppt->gc_max = GC_MAX_DEFAULT;

  ppt->degenerate = 0;          /* FIXME: may be converted to
                                   assay_type */
  ppt->detect_two_cytosines = 0;        /* 0: select the best probeset
                                           from the forward and
                                           reverse strands ; 1 two and
                                           print into a single line;
                                           2: two and print into two
                                           separate lines.  The
                                           default is to select a
                                           single best probeset for
                                           each CpG dimer */

  ppt->klass = INFINIUM_METHYL_BISULFITE;
  ppt->probe_type = INF_ASPE;
  ppt->s_strand = CODING;       /* It must be CODING.  However, no
                                   validation will be done.  The
                                   input file is responsible for the
                                   correctness of this piece of
                                   information. */

  ppt->cpg_coord = -1;
  ppt->cpg_pos = +1;

  ppt->alleles[0] = 'C';
  ppt->alleles[1] = 'G';
  ppt->alleles[2] = '\0';
  ppt->alleles[3] = '\0';
  ppt->alleles[4] = '\0';

  ppt->allele_strand = FORWARD; /* Will be reset as necessary */
  ppt->allele_coord = -1;

  ppt->allele = 0;
  ppt->allele_type = UNKNOWN;
  ppt->tss_coord = 0;

  ppt->seq = NULL;
  ppt->ori_seq = NULL;
  ppt->seq_c = NULL;
  ppt->seq_t = NULL;
  ppt->seq_rc_c = NULL;
  ppt->seq_rc_t = NULL;

  ppt->min_score_to_print = MIN_SCORE_TO_PRINT_DEFAULT;

  ppt->final_validation = 1;    /* Validation */

  return ppt;
}

/* The seqs in ppt are freed elsewhere - this makes ppt reusable */
void
probe_param_free (probe_param_t * ppt)
{
  if (ppt)
    free (ppt);
}

int
probe_param_validated (probe_param_t * ppt, const int print_p)
{
  int is_verbose;
  time_t t;
  struct utsname un;
  char hname[100];

  if (ppt->len_min < LEN_MIN_DEFAULT || ppt->len_max > LEN_MAX_DEFAULT)
    {
      fprintf (stderr, "Sequence length limit failed validation.\n");
      return 0;
    }

  assert (ppt->allele_strand == BOTH || ppt->allele_strand == FORWARD
          || ppt->allele_strand == REVERSE);

  assert (ppt->detect_two_cytosines >= 0 && ppt->detect_two_cytosines <= 2);

  if (print_p == 1)
    {
      is_verbose = get_verbose ();
      set_verbose (1);

      gethostname (hname, 100);
      time (&t);
      uname (&un);

      info0 (__FILE__, __FUNCTION__, __LINE__, "Parameters");

      info_s ("Program Package", PACKAGE);
      info_s ("Program Version", VERSION);
      info_s ("Program Author", PACKAGE_BUGREPORT);
      info_s ("Machine Architecture", un.machine);
      info_s ("Compile Machine", hname);
      info_s ("Compile Time", ctime (&t));

      info_i ("Minimum Probe Length", ppt->len_min);
      info_i ("Maximum Probe Length", ppt->len_max);
      info_i ("Minimum GC Content", ppt->gc_min);
      info_i ("Maximum GC Content", ppt->gc_max);
      info_i ("Degenerate", ppt->degenerate);
      info_i ("Detect Two Cs", ppt->detect_two_cytosines);
      info_f ("Score Threadhold", ppt->min_score_to_print);
      info_i ("Final Validation", ppt->final_validation);
#ifdef DEBUG
      info_s ("Debug", "On");
#else
      info_s ("Debug", "Off");
#endif
#ifdef CAUTION
      info_s ("Assertion", "On");
#else
      info_s ("Assertion", "Off");
#endif

      set_verbose (is_verbose);
    }

  return 1;
}

/*
  Local Variables:
  c-file-style: "gnu"
  End:
 */
