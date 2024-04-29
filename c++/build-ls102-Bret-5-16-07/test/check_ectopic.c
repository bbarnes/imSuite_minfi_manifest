
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include <check.h>

#include "check_all.h"
#include "ectopic.h"

#define TEST_FILE "../../test/data/50mer-Ectopic.txt"

static void
setup (void)
{
}

static void
teardown (void)
{
  ;
}

START_TEST (test_ectopic_count)
{
  FILE *fp;
  char seq[51];
  int ectopic, n, i;

  fp = fopen (TEST_FILE, "r");
  if (fp == NULL)
    {
      fprintf (stderr, "Ectopic test data file %s does not exist\n",
               TEST_FILE);
      exit (EXIT_FAILURE);
    }

  i = 0;
  while ((n = fscanf (fp, "%s\t%d ", seq, &ectopic)) != EOF)
    {
      if (n != 2)
        {
          fprintf (stderr, "Corrupted test data file!\n");
          continue;
        }
      fail_unless (lpm (seq, 50) == ectopic, NULL);
      i++;
    }
#ifdef DEBUG
  fprintf (stderr, "Number of ectopic tested entries: %d - all passed!\n", i);
#endif
  fclose (fp);
}
END_TEST

Suite * make_ectopic_suite (void)
{
  Suite *s;
  TCase *tc_core;

  s = suite_create ("Ectopic");
  tc_core = tcase_create ("Core Tests");
  suite_add_tcase (s, tc_core);
  tcase_add_checked_fixture (tc_core, setup, teardown);

  tcase_add_test (tc_core, test_ectopic_count);

  return s;
}

/*
  Local Variables:
  c-file-style: "gnu"
  End:
 */
