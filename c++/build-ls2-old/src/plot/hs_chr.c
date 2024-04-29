/* A small utility for enumerating Homo sapiens chromosomes */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "hs_chr.h"

int
hs_chr_to_num (const char *chr)
{
  int ret;
  ret = -1;
  if (chr == NULL || chr[0] == '\0')
    return ret;
  if (strcasecmp ("X", chr) == 0)
    ret = CHR_X;
  else if (strcasecmp ("Y", chr) == 0)
    ret = CHR_Y;
  else if (strcasecmp ("MT", chr) == 0)
    ret = CHR_MT;
  else if (strcasecmp ("UN", chr) == 0)
    ret = CHR_UN;
  else if (isdigit (chr[0]))
    {
      ret = atoi (chr);
      if (ret < 1 || ret > 22)
        ret = -1;
    }
  return ret;
}

char *
hs_num_to_chr (const int chr)
{
  static char chr_s[3];         /* ! */
  if (chr <= 22)
    sprintf (chr_s, "%d", chr);
  else if (chr == CHR_MT)
    sprintf (chr_s, "Mt");
  else if (chr == CHR_X)
    sprintf (chr_s, "X");
  else if (chr == CHR_Y)
    sprintf (chr_s, "Y");
  else if (chr == CHR_UN)
    sprintf (chr_s, "Un");
  else
    chr_s[0] = '\0';
  return chr_s;
}

/*
  Local Variables:
  c-file-style: "gnu"
  End:
 */
