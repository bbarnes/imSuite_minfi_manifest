
/* by Lixin Zhou */

#include <assert.h>

#if HAVE_CONFIG_H
#include "config.h"
#endif

#include "biolib.h"

#include "imparam.h"
#include "improbe.h"

#include "imscore-cpg.h"
#include "imscore-dup.h"
#include "imscore-gc.h"
#include "imscore-length.h"
#include "imscore-run.h"
#include "imscore-self.h"
#include "imscore-tm.h"
#include "imscore-tango.h"

/* These should be in consistence with the thermod.c implementation */
#define CONC_DNA 50.0
#define CONC_SALT 50.0

/* Perhaps we should NOT use nearest neighbor algorithm computing Tm
   since the length of the probes is large?  Considering the highly
   biased composition in methylation probes, NN may be a better
   algorithm. */
void
probe_properties (probe_t * tmp)
{
  tmp->nmer = nmer_count (tmp->probe, tmp->length);
  tmp->tm = tm_long (tmp->probe, CONC_SALT);    /* = tm_nb (tmp->probe, CONC_DNA, CONC_SALT);  */
  tmp->gc_percent = content_gc (tmp->probe) * 100.0;
  tmp->self_any = probe_self_any (tmp->probe, tmp->length);
  tmp->n_run = longest_run (tmp->probe);
  tmp->n_cpg = count_cpg (tmp->seq_ori);
  tmp->tango_count = tango_count (tmp->probe, tmp->length);

  /* Correction must be performed, now */
  if (tmp->allele_type == ALLELE_C && tmp->n_cpg > 0)
    (tmp->n_cpg)--;

#ifdef DEBUG
  debug0 (__FILE__, __FUNCTION__, __LINE__, "Probe Properties");
  debug_i ("Length", tmp->length);
  debug_i ("Nmer", tmp->nmer);
  debug_f ("Tm", tmp->tm);
  debug_f ("GC-%", tmp->gc_percent);
  debug_i ("Self-Any", tmp->self_any);
  debug_i ("N Run", tmp->n_run);
  debug_i ("N CpG", tmp->n_cpg);
  debug_i ("N Tango", tmp->tango_count);
#endif
}

void
probe_scores (probe_t * tmp)
{
  tmp->length_score = probe_length_score (tmp->length);
  tmp->nmer_score = nmer_score (tmp->nmer, tmp->allele_type);
  tmp->tm_score = tm_score (tmp->tm);
  tmp->gc_score = gc_score (tmp->gc_percent, tmp->allele_type);
  tmp->self_any_score = probe_self_any_score (tmp->self_any);
  tmp->n_run_score = longest_run_score (tmp->n_run);
  tmp->cpg_score = underlying_cpg_score (tmp->n_cpg, tmp->allele_type);
  tmp->tango_score = tango_score (tmp->tango_count);

#ifdef DEBUG
  debug0 (__FILE__, __FUNCTION__, __LINE__, "Probe Scores");
  debug_f ("Length Score", tmp->length_score);
  debug_f ("Nmer Score", tmp->nmer_score);
  debug_f ("Tm Score", tmp->tm_score);
  debug_f ("GC Score", tmp->gc_score);
  debug_f ("Self-Any Score", tmp->self_any_score);
  debug_f ("N-Run Score", tmp->n_run_score);
  debug_f ("CpG Score", tmp->cpg_score);
  debug_f ("Tango Score", tmp->tango_score);
#endif

  /* Final score if the product of all the individual scores */
  tmp->score = tmp->length_score * tmp->nmer_score * tmp->tm_score *
    tmp->gc_score * tmp->self_any_score * tmp->n_run_score * tmp->cpg_score *
    tmp->tango_score;

#ifdef DEBUG
  debug_f ("FINAL Score", tmp->score);
#endif

  assert (tmp->score >= 0.0 && tmp->score <= 1.0);
}

/*
  Local Variables:
  c-file-style: "gnu"
  End:
*/
