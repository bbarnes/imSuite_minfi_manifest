
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include <check.h>

#include "check_all.h"
#include "seqcomp.h"
#include "imscore-cpg.h"
#include "imdefs.h"
#include "cpg-utils.h"

#define TEST_FILE "../../test/data/Strand-Seq-CpG_count-CpG_min_dist.txt"

static void
setup (void)
{
}

static void
teardown (void)
{
  ;
}

START_TEST (test_cpg_count_and_mindist)
{
  FILE *fp;
  char seq[51], ori;
  int cnt, dist, n, test_cnt, test_dist, i;

  fp = fopen (TEST_FILE, "r");
  if (fp == NULL)
    {
      fprintf (stderr, "CpG test data file %s does not exist\n", TEST_FILE);
      exit (EXIT_FAILURE);
    }

  i = 0;
  while ((n = fscanf (fp, "%c\t%s\t%d\t%d ", &ori, seq, &cnt, &dist)) != EOF)
    {
      if (n != 4)
        {
          fprintf (stderr, "Corrupted test data file!\n");
          continue;
        }

      if (ori == 'R' || ori == 'F')
        ;
      else
        {
          fprintf (stderr, "Bad Strand [col 1]: %c in the data file\n", ori);
          continue;
        }

      test_cnt = count_cpg (seq);
      fail_unless (test_cnt == cnt, NULL);

      if (ori == 'F')
        // test_dist = compute_min_cpg_dist (seq, STRAND_FORWARD);
        test_dist = compute_min_cpg_dist (seq, STRAND_CONVERTED,
                                          CPG_TOP_STRAND, ALLELE_C);
      else
        // test_dist = compute_min_cpg_dist (seq, STRAND_REVERSE);
        test_dist = compute_min_cpg_dist (seq, STRAND_CONVERTED,
                                          CPG_BOT_STRAND, ALLELE_C);
      fail_unless (test_dist = dist);

      i++;
    }
#ifdef DEBUG
  fprintf (stderr, "Number of CpG tested entries: %d - all passed!\n", i);
#endif
  fclose (fp);
}
END_TEST Suite * make_cpg_suite (void)
{
  Suite *s;
  TCase *tc_core;

  s = suite_create ("Underlying_CpG");
  tc_core = tcase_create ("Core Tests");
  suite_add_tcase (s, tc_core);
  tcase_add_checked_fixture (tc_core, setup, teardown);

  tcase_add_test (tc_core, test_cpg_count_and_mindist);

  return s;
}

/*
  Local Variables:
  c-file-style: "gnu"
  End:
 */
