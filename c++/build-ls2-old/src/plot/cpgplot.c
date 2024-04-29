
/* To potential user of this program:

   The input files, particularly the TSS file, is more than just lines
   of data that are in the required format!  You are assumed to know
   what you are doing.  When in doubt, check out the code, which is
   the best documentation I can give to you.  Perplexed?  Please
   contact <lzhou@illumina.com>. */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <libgen.h>

#include "xio.h"
#include "bioseq.h"
#include "fasta.h"

#include "hs_chr.h"
#include "tss.h"
#include "cpg-score.h"
#include "cpg-dimmer.h"
#include "cpg-islands.h"
#include "plotdriver.h"

#if HAVE_CONFIG_H
#include "config.h"
#endif

/* NCBI's Human RefSeq 36.1 */
static char *chr_lookup[25][2] = {
  {"gi|89161185|ref|NC_000001.9|NC_000001", "1"},
  {"gi|89161199|ref|NC_000002.10|NC_000002", "2"},
  {"gi|89161205|ref|NC_000003.10|NC_000003", "3"},
  {"gi|89161207|ref|NC_000004.10|NC_000004", "4"},
  {"gi|51511721|ref|NC_000005.8|NC_000005", "5"},
  {"gi|89161210|ref|NC_000006.10|NC_000006", "6"},
  {"gi|89161213|ref|NC_000007.12|NC_000007", "7"},
  {"gi|51511724|ref|NC_000008.9|NC_000008", "8"},
  {"gi|89161216|ref|NC_000009.10|NC_000009", "9"},
  {"gi|89161187|ref|NC_000010.9|NC_000010", "10"},
  {"gi|51511727|ref|NC_000011.8|NC_000011", "11"},
  {"gi|89161190|ref|NC_000012.10|NC_000012", "12"},
  {"gi|51511729|ref|NC_000013.9|NC_000013", "13"},
  {"gi|51511730|ref|NC_000014.7|NC_000014", "14"},
  {"gi|51511731|ref|NC_000015.8|NC_000015", "15"},
  {"gi|51511732|ref|NC_000016.8|NC_000016", "16"},
  {"gi|51511734|ref|NC_000017.9|NC_000017", "17"},
  {"gi|51511735|ref|NC_000018.8|NC_000018", "18"},
  {"gi|42406306|ref|NC_000019.8|NC_000019", "19"},
  {"gi|51511747|ref|NC_000020.9|NC_000020", "20"},
  {"gi|51511750|ref|NC_000021.7|NC_000021", "21"},
  {"gi|89161203|ref|NC_000022.9|NC_000022", "22"},
  {"gi|89161218|ref|NC_000023.9|NC_000023", "X"},
  {"gi|89161220|ref|NC_000024.8|NC_000024", "Y"},
  {NULL, NULL}
};

void update_bioseq_name (Bioseq * seq);

static char *program_name;

static struct option const longopts[] = {
  {"CpG-Relaxed", required_argument, NULL, 'r'},
  {"CpG-Strict", required_argument, NULL, 's'},
  {"Inf-Score", required_argument, NULL, 'i'},
  {"GG-Score", required_argument, NULL, 'g'},
  {"TSS", required_argument, NULL, 't'},
  {"verbose", no_argument, NULL, 'V'},
  {"help", no_argument, NULL, 'h'},
  {"version", no_argument, NULL, 'v'},
  {NULL, 0, NULL, 0}
};

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
    If no FILE is specified, the program takes sequences\n\
    from the stdin.\n\
\n\
  -r, --CpG-Relaxed=FILE  Relaxed CpG Islands file\n\
  -s, --CpG-Strict=FILE   Strict CpG Islands file\n\
  -i, --Inf-Score=FILE    Infinium Score file\n\
  -g, --GG-Score=FILE     GoldenGate Score file\n\
  -t, --TSS=FILE          TSS file\n\
\n\
      --help              Display this help and exit.\n\
      --version           Output version information and exit.\n\
", program_name);
      puts ("\nReport bugs to <" PACKAGE_BUGREPORT ">.");
    }
  exit (status == 0 ? EXIT_SUCCESS : EXIT_FAILURE);
}

int
main (int argc, char **argv)
{
  cpg_islands *islands_s = NULL;
  cpg_islands *islands_r = NULL;
  cpg_scores *inf_scores = NULL;
  cpg_scores *gg_scores = NULL;
  tsss *ref_tsss = NULL;
  char *fn_cpg_strict = NULL, *fn_cpg_relaxed = NULL, *fn_tss = NULL;
  char *fn_gg_score = NULL, *fn_inf_score = NULL;
  FILE *fp;
  Bioseq *seq;
  int optc, i;

  program_name = basename (argv[0]);
  while ((optc = getopt_long (argc, argv, "g:i:r:s:t:hv", longopts, NULL))
         != -1)
    switch (optc)
      {
      case 'g':
        fn_gg_score = optarg;
        break;

      case 'i':
        fn_inf_score = optarg;
        break;

      case 'r':
        fn_cpg_relaxed = optarg;
        break;

      case 's':
        fn_cpg_strict = optarg;
        break;

      case 't':
        fn_tss = optarg;
        break;

      case 'h':
        usage (EXIT_SUCCESS);
        break;

      case 'v':
        fprintf (stderr, "%s %s: %s\n", PACKAGE, VERSION, program_name);
        exit (EXIT_SUCCESS);

      default:
        usage (EXIT_SUCCESS);
      }

  if (fn_cpg_strict == NULL || fn_cpg_relaxed == NULL || fn_tss == NULL ||
      fn_gg_score == NULL || fn_inf_score == NULL)
    {
      fprintf (stderr, "Five Data Files must be all specified.\n");
      usage (EXIT_FAILURE);
    }

  if (argc == 1 || optind == argc)
    argv[argc++] = "-";

  /* Read in CpG Islands strict */
  islands_s = read_cpg_islands (fn_cpg_strict);

  /* Read in CpG Islands relaxed */
  islands_r = read_cpg_islands (fn_cpg_relaxed);

  /* Read in GG scores */
  gg_scores = read_cpg_scores (fn_gg_score);

  /* Read in Inf scores */
  inf_scores = read_cpg_scores (fn_inf_score);

  /* Read in TSS */
  ref_tsss = read_tsss (fn_tss);

  plot_header ();
  for (i = optind; i < argc; i++)
    {
      if (strcmp (argv[i], "-") == 0)
        fp = stdin;
      else
        fp = xfopen (argv[i], "r");
      while ((seq = fasta_read (fp)) != NULL)
        {
          update_bioseq_name (seq);
          plot_cpgdimmer (islands_s, islands_r, gg_scores, inf_scores,
                          ref_tsss, seq);
          bioseq_free (seq);
        }
      close (fp);
    }

  /*
     cpg_scores_free (gg_scores);
     cpg_scores_free (inf_scores);
     tsss_free (ref_tsss);
     cpg_islands_free (islands_s);
     cpg_islands_free (islands_r);
   */

  return EXIT_SUCCESS;
}

void
update_bioseq_name (Bioseq * seq)
{
  int i;
  for (i = 0; chr_lookup[i][0]; i++)
    if (strcasecmp (seq->name, chr_lookup[i][0]) == 0)
      {
        strcpy (seq->name, chr_lookup[i][1]);
        break;
      }
}

/*
  Local Variables:
  c-file-style: "gnu"
  End:
 */
