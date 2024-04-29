
#define _GNU_SOURCE

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#if HAVE_CONFIG_H
#include "config.h"
#endif

#include "xmem.h"
#include "xio.h"
#include "hs_chr.h"
#include "cpg-islands.h"

static int f_cmp_cpg_island (const void *a, const void *b);

/* Expected fields from a TSV file:
   1. id[75]
   2. chromosome
   3. start,
   4. end */
cpg_islands *
read_cpg_islands (const char *fn)
{
  FILE *fp;
  cpg_islands *ret;
  cpg_island *tmp_p;
  int i, chr;
  char *line, **toks;
  int n, buf_size, N_TOKENS = 4;

  ret = cpg_islands_new (400000);       /* Educated guess for the array size */

  fprintf (stderr, "Reading CpG_Islands...\n");
  fp = xfopen (fn, "r");
  line = NULL;
  buf_size = 0;
  toks = xmalloc (N_TOKENS * sizeof (char *));
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
      chr = hs_chr_to_num (toks[1]);
      if (chr == -1)
        {
          fprintf (stderr, "Skip line of unknown chr: %s\n", toks[1]);
          continue;
        }
      tmp_p = xmalloc (sizeof (cpg_island));
      strncpy (tmp_p->id, toks[0], sizeof (tmp_p->id));
      tmp_p->chr = chr;
      tmp_p->start = atoi (toks[2]);
      tmp_p->end = atoi (toks[3]);

      /* The minimum length of CpG islands is 200-bp */
      assert (tmp_p->start >= 1);
      assert ((tmp_p->start + 199) <= tmp_p->end);

      /* Check if there is a need to realloc, and do so if necessary */
      if (ret->size[tmp_p->chr][0] >= ret->size[tmp_p->chr][1])
        {
          ret->size[tmp_p->chr][1] *= 2;
          ret->data[tmp_p->chr] = xrealloc (ret->data[tmp_p->chr],
                                            sizeof (cpg_island *) *
                                            ret->size[tmp_p->chr][1]);
        }

      ret->data[tmp_p->chr][ret->size[tmp_p->chr][0]] = tmp_p;
      ret->size[tmp_p->chr][0] += 1;
    }
  fclose (fp);
  xfree (line);
  xfree (toks);

  for (i = 0; i < CHR_ALL; i++)
    {
      ret->data[i] = realloc (ret->data[i],
                              sizeof (cpg_island *) * ret->size[i][0]);
      if ((ret->size[i][0] > 0 && ret->data[i] == NULL)
          || (ret->size[i][0] <= 0 && ret->data[i] != NULL))
        {
          fprintf (stderr, "fail to realloc\n");
          exit (EXIT_FAILURE);
        }
      fprintf (stderr, "Chr:%s CpG Islands: %d\n",
               hs_num_to_chr (i), ret->size[i][0]);
      qsort (ret->data[i], ret->size[i][0], sizeof (cpg_island *),
             f_cmp_cpg_island);
    }

  return ret;
}

cpg_islands *
cpg_islands_new (const int alloc_step)
{
  cpg_islands *ret;
  int i;

  ret = xmalloc (sizeof (cpg_islands));
  ret->data = xmalloc (CHR_ALL * sizeof (cpg_island **));
  for (i = 0; i < CHR_ALL; i++)
    {
      ret->size[i][0] = 0;
      ret->size[i][1] = alloc_step;
      ret->data[i] = xmalloc (ret->size[i][1] * sizeof (cpg_island *));
      memset (ret->data[i], '\0', ret->size[i][1] * sizeof (cpg_island *));
    }

  return ret;
}

void
cpg_islands_free (cpg_islands * p)
{
  int i, j;
  for (i = 0; i < CHR_ALL; i++)
    {
      for (j = 0; j < p->size[i][0]; j++)
        cpg_island_free (p->data[i][j]);
      free (p->data[i]);
    }
  free (p->data);
  free (p);
}

void
cpg_island_free (cpg_island * p)
{
  xfree (p);
}

/* We sort by the start coordinates; a valid set of data should not
   have overlapping CpG islands. */
static int
f_cmp_cpg_island (const void *a, const void *b)
{
  const cpg_island **aa = (const cpg_island **) a;
  const cpg_island **bb = (const cpg_island **) b;
  return (*aa)->start - (*bb)->start;
}

/*
  Local Variables:
  c-file-style: "gnu"
  End:
 */
