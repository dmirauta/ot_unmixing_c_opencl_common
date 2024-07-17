#ifndef EXTERNAL_CONCAT
#include "shared_defines.h"
#endif

void DDMP(Signature_t sa, Signature_t sb, TraceData_t td) {
  int i = 0;
  int j = 0;
  int k = 0;

  FPN ab = sa.w[0];
  FPN bb = sb.w[0];

  for (int i = 0; i < sa.n + sb.n - 1; i++) {
    td.path_a[i] = -1;
    td.path_b[i] = -1;
    td.transfer[i] = FZERO;
  }

  while (1) {

    td.path_a[k] = i;
    td.path_b[k] = j;

    if (ab < bb) {
      bb -= ab;
      td.transfer[k] = ab;
      if (i < sa.n - 1) {
        i++;
        ab = sa.w[i];
      } else {
        break;
      }
    } else {
      ab -= bb;
      td.transfer[k] = bb;
      if (j < sb.n - 1) {
        j++;
        bb = sb.w[j];
      } else {
        break;
      }
    }

    k++;
  }
}

FPN W22(Signature_t sa, Signature_t sb, TraceData_t td) {

  DDMP(sa, sb, td);

  FPN res = FZERO;
  FPN dx = FZERO;

  for (int k = 0; k < sa.n + sb.n - 1; k++) {
    if (td.path_a[k] >= 0 && td.path_b[k] >= 0) {
      dx = sa.x[td.path_a[k]] - sb.x[td.path_b[k]];
      res += td.transfer[k] * dx * dx;
    } else {
      break;
    }
  }

  return res;
}
