
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
#include "tss.h"

static int f_cmp_tss (const void *a, const void *b);

/* Expected fields from a TSV file:
   1. ID,  a string of at most 74 characters long
   2. chromosome
   3. TSS
   4. strand ('+' or '-')

   No warnings are issued on the following undesired situations:
   (1) ID length is longer than 74 chars (to be truncated to 74)
   (2) Duplicated entries */
tsss *
read_tsss (const char *fn)
{
  FILE *fp;
  tsss *ret;
  tss *tmp_p;
  int i, chr;
  char *line, **toks;
  int n, buf_size, N_TOKENS = 4;

  ret = tsss_new (40000);       /* Educated guess for the array size */

  fprintf (stderr, "Reading TSSs...\n");
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
          fprintf (stderr, "At least N tokens expected: %d.  Ingore : %s.\n",
                   N_TOKENS, line);
          continue;
        }

      chr = hs_chr_to_num (toks[1]);
      if (chr == -1)
        {
          fprintf (stderr, "Skip line of unknown chr: %s\n", toks[1]);
          continue;
        }

      tmp_p = xmalloc (sizeof (tss));
      strncpy (tmp_p->id, toks[0], sizeof (tmp_p->id));
      tmp_p->chr = chr;
      tmp_p->coord = atoi (toks[2]);
      tmp_p->strand = toks[3][0];

      /* Note the hardcoded length for promoter region and 1st exon
         region.  Also note that one coordinate may be out of the
         3'-end of the sequence. */
      if (tmp_p->strand == '+')
        {
          tmp_p->prom = tmp_p->coord - 1500;
          tmp_p->exon = tmp_p->coord + 499;

          if (tmp_p->prom < 1)
            tmp_p->prom = 1;
        }
      else if (tmp_p->strand == '-')
        {
          tmp_p->prom = tmp_p->coord + 1500;
          tmp_p->exon = tmp_p->coord - 499;
        }
      else
        {
          fprintf (stderr, "Skip line of unknown strand: %c\n",
                   tmp_p->strand);
          tss_free (tmp_p);
          continue;
        }

      if (ret->size[tmp_p->chr][0] >= ret->size[tmp_p->chr][1])
        {
          ret->size[tmp_p->chr][1] *= 2;
          ret->data[tmp_p->chr] = xrealloc (ret->data[tmp_p->chr],
                                            sizeof (tss *) *
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
      ret->data[i] = realloc (ret->data[i], sizeof (tss *) * ret->size[i][0]);
      if ((ret->size[i][0] > 0 && ret->data[i] == NULL)
          || (ret->size[i][0] <= 0 && ret->data[i] != NULL))
        {
          fprintf (stderr, "fail to realloc\n");
          exit (EXIT_FAILURE);
        }
      fprintf (stderr, "Chr:%s TSSs: %d\n",
               hs_num_to_chr (i), ret->size[i][0]);
      qsort (ret->data[i], ret->size[i][0], sizeof (tss *), f_cmp_tss);
    }

  return ret;
}

tsss *
tsss_new (const int alloc_step)
{
  tsss *ret;
  int i;

  ret = xmalloc (sizeof (tsss));
  ret->data = xmalloc (CHR_ALL * sizeof (tss **));
  for (i = 0; i < CHR_ALL; i++)
    {
      ret->size[i][0] = 0;
      ret->size[i][1] = alloc_step;
      ret->data[i] = xmalloc (ret->size[i][1] * sizeof (tss *));
      memset (ret->data[i], '\0', ret->size[i][1] * sizeof (tss *));
    }

  return ret;
}

void
tsss_free (tsss * p)
{
  int i, j;
  for (i = 0; i < CHR_ALL; i++)
    {
      for (j = 0; j < p->size[i][0]; j++)
        tss_free (p->data[i][j]);
      free (p->data[i]);
    }
  free (p->data);
  free (p);
}

void
tss_free (tss * p)
{
  xfree (p);
}

static int
f_cmp_tss (const void *a, const void *b)
{
  const tss **aa = (const tss **) a;
  const tss **bb = (const tss **) b;
  return (*aa)->coord - (*bb)->coord;
}

/*
  Local Variables:
  c-file-style: "gnu"
  End:
 */
