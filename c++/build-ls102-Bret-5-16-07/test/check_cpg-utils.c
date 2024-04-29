
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include <check.h>

#include "check_all.h"
#include "cpg-utils.h"

#define TEST_FILE "../../test/data/Strand-Seq-CpG_count-CpG_min_dist.txt"

static void cpg_string (const char *l5, const char *l3, char *str, int *len);

static void bracketed_cpg_string (const char *l5, const char *l3, char *str,
                                  int *len);

struct cpg_seq
{
  char l5[90];
  char l3[90];
  cpg_strand ori;
};

struct cpg_seq seqs[] = {
  {
   "TACTCGGCTCT", "GCATTCGATCGT", CPG_TOP_STRAND},
  {
   "CGTTAGCATAT", "TTTATCGATCGT", CPG_BOT_STRAND},
  {
   "AGCCTACTCGG", "GCATCAGATCGA", CPG_BOT_STRAND},
  {
   "TACTGCTCTAT", "GCCCATCAGATC", CPG_TOP_STRAND},
  {
   "C", "GCCCATCAGATC", CPG_UNK_STRAND},
  {
    "GCCCATCAGATC", "C", CPG_UNK_STRAND},
  {
    "GCCCATCAGATC", "G", CPG_UNK_STRAND},
  {
    "GCCCATCAGATC", "", CPG_UNK_STRAND},
  {
   "CC", "GCCCATCAGATC", CPG_UNK_STRAND},
  {
   "GG", "GCCCATCAGATC", CPG_UNK_STRAND},
  {
   "CG", "GCCCATCAGATC", CPG_UNK_STRAND},
  {
   "GC", "GCCCATCAGATC", CPG_UNK_STRAND},
  {
   "TACTGCTCTCG", "GCAG", CPG_UNK_STRAND},
  {
   "TACTGCTCTCG]", "GCAG", CPG_UNK_STRAND},
  {
   "TACTGCTCTCG[", "GCAG", CPG_UNK_STRAND},
  {
   "TACTGCTCTC]G", "GCAG", CPG_UNK_STRAND},
  {
   "TACTGCTCTC[G", "GCAG", CPG_UNK_STRAND},
  {
   "TACTGCTCTC.G", "GCAG", CPG_UNK_STRAND},
  {
   "TACTGCTCTCG.", "GCAG", CPG_UNK_STRAND},
  {
   "TATCACTTGTAGACCTTGYGCCTGG", "TACAGAATTGGGGATACAGAAGTTTGAT",
   CPG_BOT_STRAND},
  {
   "ACTATTGGGGCCAATCCCCAACCT", "TAAACTTCACAGATTTATTTTTCCATCAA",
   CPG_BOT_STRAND},
  {
   "TTTTCGATATATGTGTGTTCTTTATAAAT", "AGGGACAGCAACCTTATGACTCATGTA",
   CPG_TOP_STRAND},
  {
   "GTTGTTCTTGAAGTTTTAAAAGCA", "CTGTGTATTTTTTAAAAGATCTT", CPG_TOP_STRAND},
  {
   "", "", CPG_UNK_STRAND}
};

static void
setup (void)
{
}

static void
teardown (void)
{
  ;
}

START_TEST (test_is_cpg_locus)
{
}

END_TEST
START_TEST (test_bracket_cpg_locus)
{
}

END_TEST
START_TEST (test_unbracket_cpg_locus)
{
}

END_TEST
START_TEST (test_mask_cpg_locus)
{
}

END_TEST
START_TEST (test_mask_cpg_locus_in_place)
{
}

END_TEST
START_TEST (test_query_cpg_strand)
{
  int i;
  char seq[200];
  int len;

  for (i = 0; strlen (seqs[i].l5) > 0; i++)
    {
      cpg_string (seqs[i].l5, seqs[i].l3, seq, &len);
      fail_unless (query_cpg_strand (seq, len) == seqs[i].ori, NULL);
    }
#ifdef DEBUG
  fprintf(stderr, "Tested CpG TOP/BOT %d times.  All passed.\n", i);
#endif
}

END_TEST
START_TEST (test_query_cpg_strand_by_bracket)
{
  char seq[200];
  int len, i;

  for (i = 0; strlen (seqs[i].l5) > 0; i++)
    {
      bracketed_cpg_string (seqs[i].l5, seqs[i].l3, seq, &len);
      fail_unless (query_cpg_strand_by_bracket (seq) == seqs[i].ori, NULL);
    }
#ifdef DEBUG
  fprintf(stderr, "Tested [CpG] TOP/BOT %d times.  All passed.\n", i);
#endif
}
END_TEST

Suite * make_cpg_utils_suite (void)
{
  Suite *s;
  TCase *tc_core;

  s = suite_create ("CpG-Utils");
  tc_core = tcase_create ("Core Tests");
  suite_add_tcase (s, tc_core);
  tcase_add_checked_fixture (tc_core, setup, teardown);

  tcase_add_test (tc_core, test_is_cpg_locus);
  tcase_add_test (tc_core, test_bracket_cpg_locus);
  tcase_add_test (tc_core, test_unbracket_cpg_locus);
  tcase_add_test (tc_core, test_mask_cpg_locus);
  tcase_add_test (tc_core, test_mask_cpg_locus_in_place);
  tcase_add_test (tc_core, test_query_cpg_strand);
  tcase_add_test (tc_core, test_query_cpg_strand_by_bracket);

  return s;
}

static void
cpg_string (const char *l5, const char *l3, char *str, int *len)
{
  *len = strlen (l5) + 1;
  sprintf (str, "%scg%s", l5, l3);
}

static void
bracketed_cpg_string (const char *l5, const char *l3, char *str, int *len)
{
  *len = strlen (l5) + 2;
  sprintf (str, "%s[CG]%s", l5, l3);
}

/*
  Local Variables:
  c-file-style: "gnu"
  End:
 */
