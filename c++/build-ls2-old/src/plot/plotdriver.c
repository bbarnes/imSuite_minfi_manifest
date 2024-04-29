
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#if HAVE_CONFIG_H
#include "config.h"
#endif

#include "bioseq.h"
#include "hs_chr.h"
#include "tss.h"
#include "cpg-score.h"
#include "cpg-dimmer.h"
#include "cpg-islands.h"
#include "cmputils.h"

#include "plotdriver.h"

static void compute_dist (const cpg_island * i, const tss * t,
                          int *p, int *x);
static int compute_overlap (const int a, const int b,
                            const int x, const int y);

void
plot_header (void)
{
  printf ("Chromosome\t" "Coordinate\t"
          "CpG-Islands-Strict\t" "CpG-Islands-Relaxed\t" "TSS\t"
          "1.5-kb-Promoter\t" "0.5-kb-First-Exon\t"
          "Overlap-btw-Promoter-Strict\t" "Overlap-btw-Exon-Strict\t"
          "Overlap-btw-Promoter-Relaxed\t" "Overlap-btw-Exon-Relaxed\t"
          "GG-FScore\t" "GG-RScore\t" "GG-MaxScore\t"
          "Inf-FScore\t" "Inf-RScore\t" "Inf-MaxScore\n");
}

/* For each CpG in seq, compute the following columns:
   - Chromosome (read from the input)
   - Coordinate (walk along the sequence for CpG)
   - CpG Islands Strict (Value = ID of the containing island)
   - CpG Islands Relaxed (Value = ID of the containing island)
   - TSS ID (or Gene ID)
   - 1.5-kb Promoter (Is CpG dimmer in 1.5-kb promoter?  Value = -dist.tss')
   - 0.5-kb Exon  (Is CpG dimmer in 1.5-kb promoter?  Value = tss.dist)
   - Overlap in bp between 1.5-kb Promoter / CpG Islands Strict
   - Overlap in bp between 0.5-kb Exon / CpG Islands Strict
   - Overlap in bp between 1.5-kb Promoter / CpG Islands Relaxed
   - Overlap in bp between 0.5-kb Exon / CpG Islands Relaxed

   - GoldenGate_F_SCORE
   - GoldenGate_R_SCORE
   - GoldenGate_MAX_SCORE
   - Infinium_F_SCORE
   - Infinium_R_SCORE
   - Infinium_MAX_SCORE
*/
void
plot_cpgdimmer (const cpg_islands * islands_s, const cpg_islands * islands_r,
                const cpg_scores * gg_scores, const cpg_scores * inf_scores,
                const tsss * ref_tsss, const Bioseq * seq)
{
  int coord, chr, j, tmp_ii;
  int start_coord = 0;
  int stop_coord = seq->length - 1;
  cpg_island **tmp_s = NULL, **tmp_r = NULL;
  tss **tmp_p = NULL, **tmp_x = NULL;
  cpg_score **tmp_g = NULL, **tmp_i = NULL;
  cpg_dimmer tmp_d;
  FILE *fp = stdout;
  int s_p, s_x, r_p, r_x;

  for (j = start_coord; j < stop_coord; j++)
    if ((*(seq->seq + j) == 'C' || *(seq->seq + j) == 'c') &&
        (*(seq->seq + j + 1) == 'G' || *(seq->seq + j + 1) == 'g'))
      {
        chr = hs_chr_to_num (seq->name);
        coord = j + 1;          /* 1-based coordinate */
        tmp_d.coord = coord;

        fprintf (fp, "%s\t%d", hs_num_to_chr (chr), tmp_d.coord);

        tmp_s = bsearch (&tmp_d, islands_s->data[chr],
                         islands_s->size[chr][0], sizeof (cpg_island *),
                         f_dimmer_vs_island);
        tmp_r = bsearch (&tmp_d, islands_r->data[chr],
                         islands_r->size[chr][0], sizeof (cpg_island *),
                         f_dimmer_vs_island);
        tmp_p = bsearch (&tmp_d, ref_tsss->data[chr],
                         ref_tsss->size[chr][0], sizeof (tss *),
                         f_dimmer_vs_promoter);
        if (tmp_p == NULL)      /* A CpG dimmer is in Promoter or Exon */
          tmp_x = bsearch (&tmp_d, ref_tsss->data[chr],
                           ref_tsss->size[chr][0], sizeof (tss *),
                           f_dimmer_vs_exon);
        else
          tmp_x = NULL;
        tmp_g = bsearch (&(tmp_d.coord), gg_scores->data[chr],
                         gg_scores->size[chr][0], sizeof (cpg_score *),
                         f_srch_coord_vs_score);
        tmp_i = bsearch (&(tmp_d.coord), inf_scores->data[chr],
                         inf_scores->size[chr][0], sizeof (cpg_score *),
                         f_srch_coord_vs_score);

        if (tmp_s && tmp_p)
          compute_dist (*tmp_s, *tmp_p, &s_p, &s_x);
        else if (tmp_s && tmp_x)
          compute_dist (*tmp_s, *tmp_x, &s_p, &s_x);
        else
          s_p = s_x = 0;

        if (tmp_r && tmp_p)
          compute_dist (*tmp_r, *tmp_p, &r_p, &r_x);
        else if (tmp_r && tmp_x)
          compute_dist (*tmp_r, *tmp_x, &r_p, &r_x);
        else
          r_p = r_x = 0;

        /* Strict CpG Islands */
        if (tmp_s)
          fprintf (fp, "\t%s.%d.%d", hs_num_to_chr (chr),
                   (*tmp_s)->start, (*tmp_s)->end);
        else
          fprintf (fp, "\t.");

        /* Relaxed CpG Islands */
        if (tmp_r)
          fprintf (fp, "\t%s.%d.%d", hs_num_to_chr (chr),
                   (*tmp_r)->start, (*tmp_r)->end);
        else
          fprintf (fp, "\t.");

        /* TSS ID */
        if (tmp_p)
          fprintf (fp, "\t%s", (*tmp_p)->id);
        else if (tmp_x)
          fprintf (fp, "\t%s", (*tmp_x)->id);
        else
          fprintf (fp, "\t.");

        /* Promoter region.  Excise caution for multiple gene entries! */
        if (tmp_p)              /* TSS = (*tmp_p)->coord */
          {
            /* coordinates: from -1500 to -1 */
            if ((*tmp_p)->strand == '+')
              fprintf (fp, "\t%d", coord - (*tmp_p)->coord);
            else if ((*tmp_p)->strand == '-')
              fprintf (fp, "\t%d", (*tmp_p)->coord - coord);
            else
              fprintf (fp, "\t?");
          }
        else
          fprintf (fp, "\t.");

        /* First exon region.  Excise caution for multiple gene entries! */
        if (tmp_x)              /* TSS = (*tmp_p)->coord */
          {
            /* coordinates: from 0 to 499 */
            if ((*tmp_x)->strand == '+')
              fprintf (fp, "\t%d", coord - (*tmp_x)->coord);
            else if ((*tmp_x)->strand == '-')
              fprintf (fp, "\t%d", (*tmp_x)->coord - coord);
            else
              fprintf (fp, "\t?");
          }
        else
          fprintf (fp, "\t.");

        /* Overlap between promter and CpG Island Strict */
        fprintf (fp, "\t%d", s_p);

        /* Overlap between exon and CpG Island Strict */
        fprintf (fp, "\t%d", s_x);

        /* Overlap between promoter and CpG Island Relax */
        fprintf (fp, "\t%d", r_p);

        /* Overlap between exon and CpG Island relax */
        fprintf (fp, "\t%d", r_x);

        /* GG Scores */
        if (tmp_g)
          fprintf (fp, "\t%.2f\t%.2f\t%.2f",
                   (*tmp_g)->f_score, (*tmp_g)->r_score, (*tmp_g)->max_score);
        else
          fprintf (fp, "\t.\t.\t.");

        /* Inf Scores */
        if (tmp_i)
          fprintf (fp, "\t%.2f\t%.2f\t%.2f",
                   (*tmp_i)->f_score, (*tmp_i)->r_score, (*tmp_i)->max_score);
        else
          fprintf (fp, "\t.\t.\t.");

        fprintf (fp, "\n");
      }
}

static void
compute_dist (const cpg_island * i, const tss * t, int *p, int *x)
{
  assert (i != NULL);
  assert (t != NULL);

  if (t->strand == '+')
    {
      *p = compute_overlap (i->start, i->end, t->prom, t->coord - 1);
      *x = compute_overlap (i->start, i->end, t->coord, t->exon);
    }
  else if (t->strand == '-')
    {
      *p = compute_overlap (i->start, i->end, t->coord + 1, t->prom);
      *x = compute_overlap (i->start, i->end, t->exon, t->coord);
    }
}

static int
compute_overlap (const int a, const int b, const int x, const int y)
{
  int ret = 0;
  int min_ax, max_by;

  assert (a <= b);
  assert (x <= y);

  if ((a < x && b < x) || (a > y && b > y))     /* No overlap */
    ret = 0;
  else
    {
      min_ax = a < x ? a : x;
      max_by = b > y ? b : y;
      ret = max_by - min_ax + 1 - abs (a - x) - abs (b - y);
    }

  return ret;
}

/*
  Local Variables:
  c-file-style: "gnu"
  End:
 */
