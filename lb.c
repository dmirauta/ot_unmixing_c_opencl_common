#ifndef EXTERNAL_CONCAT
#include "dists.h"
#include "shared_defines.h"
#endif

/// Computes a barycenter between discrete measures according to displacement
/// interpolation and barycentric projection. Expects signatures of the same
/// length (does not check). xba: positions of resulting barycenter (weights
/// will be the same as in signature a)
void left_bary(Signature_t sa, Signature_t sb, FPN *xba, FPN λ, TraceData_t td,
               FPN eps)

{
  DDMP(sa, sb, td);

  for (int i = 0; i < sa.n; i++) {
    xba[i] = λ * sa.x[i];
  }

  int i, j;
  for (int k = 0; k < sa.n * 2 - 1; k++) {
    i = td.path_a[k];
    j = td.path_b[k];
    if (i >= 0 && j >= 0) {
      xba[i] += (1 - λ) * td.transfer[k] * sb.x[j] / (sa.w[i] + eps);
    }
  }
}
