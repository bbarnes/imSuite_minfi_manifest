
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <getopt.h>
#include <libgen.h>

#include "xio.h"
#include "xmem.h"

#if HAVE_CONFIG_H
#include "config.h"
#endif

static char *program_name;

static struct option const longopts[] = {
  {"help", no_argument, NULL, 'h'},
  {"version", no_argument, NULL, 'v'},
  {NULL, 0, NULL, 0}
};

static void update_bioseq_name (char *id);

static void
usage (int status)
{
  if (status != 0)
    fprintf (stderr, "Try `%s --help' for more information.\n", program_name);
  else
    {
      printf ("\
Usage: %s [OPTION] [FILE [FILE [FILE ...]]]\n\
\n\
      --help              Display this help and exit.\n\
      --version           Output version information and exit.\n\
", program_name);
      puts ("\nReport bugs to <" PACKAGE_BUGREPORT ">.");
    }
  exit (status == 0 ? EXIT_SUCCESS : EXIT_FAILURE);
}

/* Human RefSeq 36.1 GI to chromosome lookup table */
static char *chr_lookup[25][2] = {
  {"89161185", "1"},
  {"89161199", "2"},
  {"89161205", "3"},
  {"89161207", "4"},
  {"51511721", "5"},
  {"89161210", "6"},
  {"89161213", "7"},
  {"51511724", "8"},
  {"89161216", "9"},
  {"89161187", "10"},
  {"51511727", "11"},
  {"89161190", "12"},
  {"51511729", "13"},
  {"51511730", "14"},
  {"51511731", "15"},
  {"51511732", "16"},
  {"51511734", "17"},
  {"51511735", "18"},
  {"42406306", "19"},
  {"51511747", "20"},
  {"51511750", "21"},
  {"89161203", "22"},
  {"89161218", "X"},
  {"89161220", "Y"},
  {NULL, NULL}
};

typedef float a[2];

int
main (int argc, char **argv)
{
  int optc;
  FILE *fp;
  int i, chr;
  char *line, **toks;
  int n, buf_size, N_TOKENS = 13;       /* Minimum Columns: 13 */

  int done_header = 0;

  char chr_id[50];
  int coord;
  int nn = 250000000;           /* ... largest human chromosome */
  a *dat;

  program_name = basename (argv[0]);
  while ((optc = getopt_long (argc, argv, "hv", longopts, NULL)) != -1)
    switch (optc)
      {
      case 'h':
        usage (EXIT_SUCCESS);
        break;

      case 'v':
        fprintf (stderr, "%s %s: %s\n", PACKAGE, VERSION, program_name);
        exit (EXIT_SUCCESS);

      default:
        usage (EXIT_SUCCESS);
      }

  if (argc == 1 || optind == argc)
    argv[argc++] = "-";

  dat = xcalloc (nn, sizeof (a));
  dat = memset (dat, 'A', nn * sizeof (a));     /* score <= 1.0 */

  /* Actually, we expect one chromosome in one file at a time. */
  /* for (i = optind; i < argc; i++) */
  for (i = optind; i <= optind; i++)
    {
      if (strcmp (argv[i], "-") == 0)
        fp = stdin;
      else
        fp = xfopen (argv[i], "r");

      line = NULL;
      buf_size = 0;
      toks = xmalloc (N_TOKENS * sizeof (char *));

      for (; (n = getline (&line, &buf_size, fp)) != -1;)
        {
          /* The first line is treated as the header */
          if (done_header == 0)
            {
              done_header = 1;
              continue;
            }

          chomp_line (line, n);
          if (is_comment_line (line))
            continue;

          /* col  1: Seq_ID
             col  3: CpG_Coord
             col  5: Allele_Strand
             col 12: Probeset_Score */
          if (parse_tokens (line, toks, N_TOKENS) == 1)
            {
              fprintf (stderr, "N tokens expected: %d.  Ingore : %s.\n",
                       N_TOKENS, line);
              continue;
            }

          strcpy (chr_id, toks[0]);     /* Oh, well... */
          coord = atoi (toks[2]);
          if (coord > nn)
            {
              fprintf (stderr, "Sorry, the coord is too large: %d\n", coord);
              exit (EXIT_FAILURE);
            }

          if (toks[4][0] == 'F')
            dat[coord][0] = atof (toks[11]);
          else if (toks[4][0] == 'R')
            dat[coord][1] = atof (toks[11]);
          else
            {
              fprintf (stderr, "Unknown strand: %c\n", toks[4][0]);
              exit (EXIT_FAILURE);
            }
        }

      fclose (fp);
      xfree (toks);
      xfree (line);
    }

  fprintf (stdout, "Chromosome\tCoordinate\tF-Score\tR-Score\tMax-Score\n");
  update_bioseq_name (chr_id);
  for (i = 0; i < nn; i++)
    if (dat[i][0] <= 1.0)
      {
        fprintf (stdout, "%s\t", chr_id);
        fprintf (stdout, "%d\t", i);
        fprintf (stdout, "%.2f\t%.2f\t", dat[i][0], dat[i][1]);
        fprintf (stdout, "%.2f\n",
                 dat[i][0] > dat[i][1] ? dat[i][0] : dat[i][1]);
      }

  return EXIT_SUCCESS;
}

static void
update_bioseq_name (char *id)
{
  int i;
  for (i = 0; chr_lookup[i][0]; i++)
    if (strcasecmp (id, chr_lookup[i][0]) == 0)
      {
        strcpy (id, chr_lookup[i][1]);
        break;
      }
}

/*
  Local Variables:
  c-file-style: "gnu"
  End:
 */
