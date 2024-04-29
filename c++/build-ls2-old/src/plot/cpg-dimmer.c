
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
#include "cpg-dimmer.h"

static int f_cmp_cpg_dimmer (const void *a, const void *b);

/* Expected fields from a TSV file:
   1. chromosome
   2. coordinate */
cpg_dimmers *
read_cpg_dimers (const char *fn)
{
  FILE *fp;
  cpg_dimmers *ret;
  cpg_dimmer *tmp_p;
  int i, chr;
  char *line, **toks;
  int n, buf_size, N_TOKENS = 2;

  ret = cpg_dimmers_new (3000000);      /* Educated guess for the array size */

  fprintf (stderr, "Reading CpG_Dimmers...\n");
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
      tmp_p = xmalloc (sizeof (cpg_dimmer));
      tmp_p->chr = chr;
      tmp_p->coord = atoi (toks[1]);
      tmp_p->strand = '+';      /* ! */

      /* Check if there is a need to realloc, and do so if necessary */
      if (ret->size[tmp_p->chr][0] >= ret->size[tmp_p->chr][1])
        {
          ret->size[tmp_p->chr][1] *= 2;
          ret->data[tmp_p->chr] = xrealloc (ret->data[tmp_p->chr],
                                            sizeof (cpg_dimmer *) *
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
                              sizeof (cpg_dimmer *) * ret->size[i][0]);
      if ((ret->size[i][0] > 0 && ret->data[i] == NULL)
          || (ret->size[i][0] <= 0 && ret->data[i] != NULL))
        {
          fprintf (stderr, "fail to realloc\n");
          exit (EXIT_FAILURE);
        }
      fprintf (stderr, "Chr:%s CpG Dimmers: %d\n",
               hs_num_to_chr (i), ret->size[i][0]);
      qsort (ret->data[i], ret->size[i][0], sizeof (cpg_dimmer *),
             f_cmp_cpg_dimmer);
    }

  return ret;
}

cpg_dimmers *
cpg_dimmers_new (const int alloc_step)
{
  cpg_dimmers *ret;
  int i;

  ret = xmalloc (sizeof (cpg_dimmers));
  ret->data = xmalloc (CHR_ALL * sizeof (cpg_dimmer **));
  for (i = 0; i < CHR_ALL; i++)
    {
      ret->size[i][0] = 0;
      ret->size[i][1] = alloc_step;
      ret->data[i] = xmalloc (ret->size[i][1] * sizeof (cpg_dimmer *));
      memset (ret->data[i], '\0', ret->size[i][1] * sizeof (cpg_dimmer *));
    }

  return ret;
}

void
cpg_dimmers_free (cpg_dimmers * p)
{
  int i, j;
  for (i = 0; i < CHR_ALL; i++)
    {
      for (j = 0; j < p->size[i][0]; j++)
        cpg_dimmer_free (p->data[i][j]);
      free (p->data[i]);
    }
  free (p->data);
  free (p);
}

void
cpg_dimmer_free (cpg_dimmer * p)
{
  xfree (p);
}

static int
f_cmp_cpg_dimmer (const void *a, const void *b)
{
  const cpg_dimmer **aa = (const cpg_dimmer **) a;
  const cpg_dimmer **bb = (const cpg_dimmer **) b;
  return (*aa)->coord - (*bb)->coord;
}

/*
  Local Variables:
  c-file-style: "gnu"
  End:
 */
