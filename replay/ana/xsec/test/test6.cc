/* nag_bessel_j1 (s17afc) Example Program.
 *
 * CLL6I261D/CLL6I261DL Version.
 *
 * Copyright 2017 Numerical Algorithms Group.
 *
 * Mark 26.1, 2017.
 */

#include <nag.h>
#include <stdio.h>
#include <nag_stdlib.h>
#include <nags.h>

int main(void)
{
  Integer exit_status = 0;
  double x, y;
  NagError fail;

  INIT_FAIL(fail);

  /* Skip heading in data file */
  scanf("%*[^\n]");
  printf("nag_bessel_j1 (s17afc) Example Program Results\n");
  printf("     x           y\n");
  while (scanf("%lf", &x) != EOF)
  {
    /* nag_bessel_j1 (s17afc).
     * Bessel function J_1(x)
     */
    y = nag_bessel_j1(x, &fail);
    if (fail.code != NE_NOERROR) {
      printf("Error from nag_bessel_j1 (s17afc).\n%s\n", fail.message);
      exit_status = 1;
      goto END;
    }
    printf("%12.3e%12.3e\n", x, y);
  }

END:
  return exit_status;
}
