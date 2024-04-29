
/* by Lixin Zhou, April 2006 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include <getopt.h>
#include <libgen.h>

#include "biolib.h"

#include "imdefs.h"
#include "imparam.h"
#include "imdriver.h"

#include "imscore-dup.h"
#include "imscore-tango.h"

#if HAVE_CONFIG_H
#include "config.h"
#endif

/* The input file should have the following minimum piece of information:
   - seq-name
   - tss-coord
   - sequence
*/

static char *program_name;

static struct option const longopts[] = {
  {"tango-file", required_argument, NULL, 'a'},
  {"13mer-file", required_argument, NULL, 'n'},
  {"fasta-format", no_argument, NULL, 'f'},
  {"asaay-format", required_argument, NULL, 'o'},
  {"two-probesets-per-CpG", required_argument, NULL, 't'},
  {"degenerate", no_argument, NULL, 'd'},
  {"score-cutoff", required_argument, NULL, 's'},
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
    from the stdin.  The format is TSV: seq-name, sequence\n\
\n\
  -f, --fasta-format      The default input file is TSV; if this switch\n\
                          is on, the input file is in Fasta format.\n\
\n\
  -o, --assay-format      The default is ASPE\n\
                          Valid: ASPE, SBE, [more to be added].\n\
\n\
  -t, --two-probesets-per-CpG=[0,1,2]\n\
                          0: The default is to select the best probeset\n\
                          1: Both probesets are printed into a single line\n\
                          2: Both probesets are printed into separate lines\n\
\n\
  -d, --degenerate        Design degenerate probes.  Default: regular\n\
\n\
  -s, --score-cutoff=s    Print probes whose final score >= \"s\"\n\
\n\
  -n, --13mer-file=FILE   13-mer library file in binary format\n\
\n\
  -a, --tango-file=FILE   11-mer Tango file in binary format\n\
\n\
  -V, --verbose           More information is printed to stderr.\n\
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
  probe_param_t *ppt;
  char **tokens, *line, *fn_13mer_b = NULL, *fn_tango_b = NULL;
  int optc, i, ntokens, is_fasta = 0, is_verbose = 0;
  FILE *fp;

  set_verbose (is_verbose);
  ppt = probe_param_new ();

  program_name = basename (argv[0]);
  while ((optc = getopt_long (argc, argv, "a:dfn:o:s:t:Vhv", longopts, NULL))
         != -1)
    switch (optc)
      {
      case 'a':
        fn_tango_b = optarg;
        break;

      case 'd':
        ppt->degenerate = 1;
        break;

      case 'f':
        is_fasta = 1;
        break;

      case 'n':
        fn_13mer_b = optarg;
        break;

      case 'o':
        if (strcasecmp (optarg, "ASPE") == 0)
          ppt->probe_type = INF_ASPE;
        else if (strcasecmp (optarg, "SBE") == 0)
          ppt->probe_type = INF_SBE;
        else
          {
            fprintf (stderr, "Not supported assay format: %s.\n", optarg);
            usage (EXIT_FAILURE);
          }
        break;

      case 's':
        ppt->min_score_to_print = (float) atof (optarg);
        break;

      case 't':
        ppt->detect_two_cytosines = atoi (optarg);
        if (ppt->detect_two_cytosines < 0 || ppt->detect_two_cytosines > 2)
          {
            fprintf (stderr, "Value of t must be 0, 1 or 2.\n");
            usage (EXIT_FAILURE);
          }
        break;

      case 'V':
        is_verbose = 1;
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

  set_verbose (is_verbose);

  if (!probe_param_validated (ppt, 1))
    {
      fprintf (stderr, "Please check the parameters and try again.\n");
      usage (EXIT_FAILURE);
    }

  if (argc == 1 || optind == argc)
    argv[argc++] = "-";

  info0 (__FILE__, __FUNCTION__, __LINE__, "Load 13mer library...");
  info_s ("Library file", fn_13mer_b);

  read_tango_11mer_data (fn_tango_b, 11);
  read_13mer_data (fn_13mer_b, 13);

  info_s ("Done", "Loading 13mer library");

  /* FIXME: The coordinate of TSS is to be provided in the input file */
  for (i = optind; i < argc; i++)
    {
      info0 (__FILE__, __FUNCTION__, __LINE__,
             "Read and process sequence files");
      info_s ("Process sequence file", argv[i]);

      if (strcmp (argv[i], "-") == 0)
        fp = stdin;
      else
        fp = xfopen (argv[i], "r");

      if (is_fasta)
        {
          while ((ppt->ori_seq = fasta_read (fp)) != NULL)
            {
              info0 (__FILE__, __FUNCTION__, __LINE__, "Sequence info");
              info_s ("Sequence Name", ppt->ori_seq->name);
              info_i ("Sequence Length", ppt->ori_seq->length);

              ppt->tss_coord = ppt->ori_seq->length - 30 - 1;
              design_probeset_seq (ppt);
              bioseq_free (ppt->ori_seq);
            }
        }

      else
        {
          for (; (line = (char *) read_a_line (fp)) && !feof (fp);)
            {
              if (!is_comment_line (line))
                {
                  tokens = tokenize (line, "\t", 0, &ntokens);
                  ppt->ori_seq = bioseq_new (tokens[0], NULL, tokens[1]);

                  info0 (__FILE__, __FUNCTION__, __LINE__, "Sequence info");
                  info_s ("Sequence Name", ppt->ori_seq->name);
                  info_i ("Sequence Length", ppt->ori_seq->length);

                  ppt->tss_coord = ppt->ori_seq->length - 30 - 1;
                  design_probeset_seq (ppt);
                  bioseq_free (ppt->ori_seq);
                }
              free (line);
              delete_tokens (tokens, ntokens);
            }
          xfree (line);
        }

      if (fp != stdin)
        fclose (fp);

      info_s ("Done", "Processing a sequence file");
    }

  probe_param_free (ppt);
  delete_13mer_data ();

  return EXIT_SUCCESS;
}

/*
  Local Variables:
  c-file-style: "gnu"
  End:
 */
