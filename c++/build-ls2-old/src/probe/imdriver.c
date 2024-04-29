
/* by Lixin Zhou */

#include <assert.h>
#include <stdlib.h>
#include <ctype.h>

#if HAVE_CONFIG_H
#include "config.h"
#endif

#include "biolib.h"

#include "imparam.h"
#include "improbeset.h"
#include "imdriver.h"

static Bioseq *bioseq_new1 (const Bioseq * in);

void
design_probeset_seq (probe_param_t * ppt)
{
  probesets_t *psets = NULL;
  char *tmp;
  long int coord, start_coord, stop_coord, i;

  info0 (__FILE__, __FUNCTION__, __LINE__,
         "Prepare the template sequences...");

  /* Covert the input sequence into uppcase */
  assert (ppt->ori_seq->seq != NULL);
  for (i = 0; ppt->ori_seq->seq[i]; i++)
    ppt->ori_seq->seq[i] = toupper (ppt->ori_seq->seq[i]);

  /* Prepare four template sequences.  It saves time by doing all
     these at one time at a cost of (hopefully small) space.  */
  ppt->seq_c = bioseq_new1 (ppt->ori_seq);
  ppt->seq_t = bioseq_new1 (ppt->ori_seq);
  ppt->seq_rc_c = bioseq_new1 (ppt->ori_seq);
  ppt->seq_rc_t = bioseq_new1 (ppt->ori_seq);

  tmp = xmalloc (ppt->ori_seq->length + 1);
  nuseq_rc (ppt->ori_seq->seq, ppt->ori_seq->length, tmp);

  if (ppt->degenerate)
    {
      bisulfite (ppt->ori_seq->seq, ppt->seq_c->seq);
      bisulfite (ppt->ori_seq->seq, ppt->seq_t->seq);
      bisulfite (tmp, ppt->seq_rc_c->seq);
      bisulfite (tmp, ppt->seq_rc_t->seq);
    }
  else
    {
      bisulfite_methyl (ppt->ori_seq->seq, ppt->seq_c->seq);
      bisulfite_unmethyl (ppt->ori_seq->seq, ppt->seq_t->seq);
      bisulfite_methyl (tmp, ppt->seq_rc_c->seq);
      bisulfite_unmethyl (tmp, ppt->seq_rc_t->seq);
    }
  xfree (tmp);

  info_s ("Done", "Preparing the template sequences.");
  info0 (__FILE__, __FUNCTION__, __LINE__, "Design probes...");

  /* CpG near the end is not simply ignored.  It's taken care of
     later, if a CpG is too close to an end of the sequence... */
  start_coord = 0 + 5;
  stop_coord = ppt->ori_seq->length - 1 - 5;

  i = 0;
  for (coord = start_coord; coord < stop_coord; coord++)
    {
      if ((*(ppt->ori_seq->seq + coord) == 'C' ||
           *(ppt->ori_seq->seq + coord) == 'c') &&
          (*(ppt->ori_seq->seq + coord + 1) == 'G' ||
           *(ppt->ori_seq->seq + coord + 1) == 'g'))
        {
          ppt->cpg_coord = coord;       /* 0-based index; '+' strand */
          ppt->allele_coord = coord;    /* Changed later for '-' strand */

          /* TSS = +1; immediately upstream: -1.  Note that both coord
             and tss_coord are string indexes, so that the relative
             position, cpg_pos, is correct */
          ppt->cpg_pos = coord - ppt->tss_coord;
          if (ppt->cpg_pos >= 0)
            (ppt->cpg_pos)++;

#ifdef DEBUG
          debug0 (__FILE__, __FUNCTION__, __LINE__,
                  "Design probes at this position:");
          debug_i ("CpG index", ppt->cpg_coord);
          debug_i ("Total CpG Processed:", ++i);
#endif

          psets = design_probesets (ppt);
          print_probesets (ppt, psets);

          probesets_free (psets);
        }
    }

  info_s ("Done", "Designing probes.");
  info0 (__FILE__, __FUNCTION__, __LINE__, "Delete template sequences...");

  bioseq_free (ppt->seq_c);
  bioseq_free (ppt->seq_t);
  bioseq_free (ppt->seq_rc_c);
  bioseq_free (ppt->seq_rc_t);

  info_s ("Done", "Deleting template sequences.");
}

static Bioseq *
bioseq_new1 (const Bioseq * in)
{
  Bioseq *out;

  assert (in != NULL);
  assert (in->seq != NULL);

  out = bioseq_new (in->name, NULL, NULL);
  out->seq = (char *) xmalloc (in->length + 1);
  out->length = in->length;

  return out;
}

/*
  Local Variables:
  c-file-style: " gnu "
  End:
 */
