
#ifndef EXTERNAL_CONCAT
#include "shared_defines.h"
#endif

FPN segment_cost(const FPN *xc, const FPN *c, GLB const FPN *XBfs,
                 GLB const FPN *WBfs, // passed through flat and sorted
                 FPN *av, FPN *cts, GLB const int *karr) {
  int i, j, k;

  // concatenate basis signatures
  for (i = 0; i < K; i++) {
    for (j = 0; j < N; j++) {
      k = j + i * N;
      cts[k] = av[karr[k]] * WBfs[k];
    }
  }

  i = 0;
  j = 0;

  // slacks
  FPN sc = c[0];
  FPN sct = cts[0];
  FPN W = FZERO;
  FPN dx = FZERO;

  int max_iter = (K + 1) * N;
  for (k = 0; k < max_iter; k++) {

    dx = xc[i] - XBfs[j];
    if (sc < sct) {
      sct -= sc;
      W += sc * dx * dx;
      if (i < N - 1) {
        i += 1;
        sc = c[i];
      } else {
        break;
      }
    } else {
      sc -= sct;
      W += sct * dx * dx;
      if (j < N * K - 1) {
        j += 1;
        sct = cts[j];
      } else {
        break;
      }
    }
  }

  return W;
}

FPN segment_cost_aug(const FPN *xc, const FPN *c, GLB const FPN *XBfs,
                     GLB const FPN *WBfs, // passed through flat and sorted
                     FPN *av, FPN *cts, GLB const int *karr) {
  FPN sum = FZERO;
  FPN sum_squares = FZERO;
  for (int i = 0; i < K; i++) {
    if (av[i] < 0) {
      return FLT_MAX;
    }
    sum += av[i];
    sum_squares += av[i] * av[i];
  }
  FPN sc = segment_cost(xc, c, XBfs, WBfs, av, cts, karr);
  sc += 10 * (sum - 1) * (sum - 1); // non-normalisation penalty
  sc -= SPARSITY_PROMOTION * sum_squares;
  return sc;
}

FPN segment_cost_rc(const FPN *xc, const FPN *c, GLB const FPN *XBfs,
                    GLB const FPN *WBfs, // passed through flat and sorted
                    FPN *av, FPN *cts, GLB const int *karr, FPN *Bt) {
  int i, j, k;

  // concatenate basis signatures
  for (i = 0; i < K; i++) {
    for (j = 0; j < N; j++) {
      k = j + i * N;
      cts[k] = av[karr[k]] * WBfs[k];
      Bt[k] = FZERO;
    }
  }

  i = 0;
  j = 0;

  // slacks
  FPN sc = c[0];
  FPN sct = cts[0];
  FPN W = FZERO;
  FPN dx = FZERO;

  int max_iter = (K + 1) * N;
  for (k = 0; k < max_iter; k++) {

    dx = xc[i] - XBfs[j];
    if (sc < sct) {
      sct -= sc;
      Bt[karr[j] * N + i] += sc;
      W += sc * dx * dx;
      if (i < N - 1) {
        i += 1;
        sc = c[i];
      } else {
        break;
      }
    } else {
      sc -= sct;
      Bt[karr[j] * N + i] += sct;
      W += sct * dx * dx;
      if (j < N * K - 1) {
        j += 1;
        sct = cts[j];
      } else {
        break;
      }
    }
  }

  return W;
}
