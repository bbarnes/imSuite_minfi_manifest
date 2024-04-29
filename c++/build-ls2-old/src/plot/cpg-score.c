
#define _GNU_SOURCE

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#if HAVE_CONFIG_H
#include "config.h"
#endif

#include "xmem.h"
#include "xio.h"
#include "hs_chr.h"
#include "cpg-score.h"

static int f_cmp_cpg_score (const void *a, const void *b);

/* Expected fields from a TSV file:
   1. chromosome
   2. coordinate
   3. f_score
   4. r_score
   5. max_score -- may be made optional
 */
cpg_scores *
read_cpg_scores (const char *fn)
{
  FILE *fp;
  cpg_scores *ret;
  cpg_score *tmp_p;
  int i, chr;
  char *line, **toks;
  int n, buf_size, N_TOKENS = 5;

  ret = cpg_scores_new (3000000);       /* Educated guess for the array size */

  fprintf (stderr, "Reading CpG_Scores...\n");
  line = NULL;
  buf_size = 0;
  toks = xmalloc (N_TOKENS * sizeof (char *));
  fp = xfopen (fn, "r");
  for (; (n = getline (&line, &buf_size, fp)) != -1;)
    {
      chomp_line (line, n);
      if (is_comment_line (line))
        continue;

      if (parse_tokens (line, toks, N_TOKENS) == 1)
        {
          fprintf (stderr, "N tokens expected: %d.  Ingore : %s.\n",
                   N_TOKENS, line);
          continue;
        }

      /* Get and validate the data */
      chr = hs_chr_to_num (toks[0]);
      if (chr == -1)
        {
          fprintf (stderr, "Skip line of unknown chr: %s\n", toks[0]);
          continue;
        }
      tmp_p = xmalloc (sizeof (cpg_score));
      tmp_p->chr = chr;
      tmp_p->coord = atoi (toks[1]);
      tmp_p->f_score = atof (toks[2]);
      tmp_p->r_score = atof (toks[3]);
      tmp_p->max_score = atof (toks[4]);

      /* Check if there is a need to realloc, and do so if necessary */
      if (ret->size[tmp_p->chr][0] >= ret->size[tmp_p->chr][1])
        {
          ret->size[tmp_p->chr][1] *= 2;
          ret->data[tmp_p->chr] = xrealloc (ret->data[tmp_p->chr],
                                            sizeof (cpg_score *) *
                                            ret->size[tmp_p->chr][1]);
        }

      ret->data[tmp_p->chr][ret->size[tmp_p->chr][0]] = tmp_p;
      ret->size[tmp_p->chr][0] += 1;
    }
  fclose (fp);
  xfree (toks);
  xfree (line);

  for (i = 0; i < CHR_ALL; i++)
    {
      ret->data[i] = realloc (ret->data[i],
                              sizeof (cpg_score *) * ret->size[i][0]);
      if ((ret->size[i][0] > 0 && ret->data[i] == NULL)
          || (ret->size[i][0] <= 0 && ret->data[i] != NULL))
        {
          fprintf (stderr, "fail to realloc\n");
          exit (EXIT_FAILURE);
        }
      fprintf (stderr, "Chr:%s CpG Scores: %d\n",
               hs_num_to_chr (i), ret->size[i][0]);
      qsort (ret->data[i], ret->size[i][0], sizeof (cpg_score *),
             f_cmp_cpg_score);
    }

  return ret;
}

cpg_scores *
cpg_scores_new (const int alloc_step)
{
  cpg_scores *ret;
  int i;

  ret = xmalloc (sizeof (cpg_scores));
  ret->data = xmalloc (CHR_ALL * sizeof (cpg_score **));
  for (i = 0; i < CHR_ALL; i++)
    {
      ret->size[i][0] = 0;
      ret->size[i][1] = alloc_step;
      ret->data[i] = xmalloc (ret->size[i][1] * sizeof (cpg_score *));
      memset (ret->data[i], '\0', ret->size[i][1] * sizeof (cpg_score *));
    }

  return ret;
}

void
cpg_scores_free (cpg_scores * p)
{
  int i, j;
  for (i = 0; i < CHR_ALL; i++)
    {
      for (j = 0; j < p->size[i][0]; j++)
        cpg_score_free (p->data[i][j]);
      free (p->data[i]);
    }
  free (p->data);
  free (p);
}

void
cpg_score_free (cpg_score * p)
{
  xfree (p);
}

int
f_srch_coord_vs_score (const void *a, const void *b)
{
  const int *coord = (const int *) a;
  const cpg_score **bb = (const cpg_score **) b;
  return (*coord) - (*bb)->coord;
}

static int
f_cmp_cpg_score (const void *a, const void *b)
{
  const cpg_score **aa = (const cpg_score **) a;
  const cpg_score **bb = (const cpg_score **) b;
  return (*aa)->coord - (*bb)->coord;
}

/*
  Local Variables:
  c-file-style: "gnu"
  End:
 */
