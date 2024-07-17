
#ifndef EXTERNAL_CONCAT
#include "quicksort.c"
#include "shared_defines.h"
#endif

// NOTE: Still currently somewhat ocl kernel specific
// TODO: Remedy above

/// maximum edge length of bounding box (of vertices in V)
/// used to measure if the simplex has shrunk sufficiently
FPN melbb(FPN V[K][K], FPN mins[K], FPN maxs[K]) {
  for (int i = 0; i < K; i++) {
    mins[i] = 2 * FONE;
    maxs[i] = -FONE;
  }

  for (int i = 0; i < K; i++) {
    for (int j = 0; j < K; j++) {
      mins[i] = fmin(mins[i], V[j][i]);
      maxs[i] = fmax(maxs[i], V[j][i]);
    }
  }

  // store diams in max vec now
  for (int i = 0; i < K; i++) {
    maxs[i] -= mins[i];
  }

  FPN diam = maxs[0];
  for (int i = 1; i < K; i++) {
    diam = fmax(diam, maxs[i]);
  }

  return diam;
}

void conv_add(FPN *a, FPN *b, FPN *c, FPN λ) {
  for (int i = 0; i < K; i++) {
    c[i] = (1 - λ) * a[i] + λ * b[i];
  }
}

void copy_vec(FPN *dst, FPN *src) {
  for (int i = 0; i < K; i++) {
    dst[i] = src[i];
  }
}

typedef struct NM_info {
  int iter; // number of iterations used
  FPN diam; // measure final simplex size for debug/sanity check purposes
} NM_info_t;

/// Optimises FPN NM_func(FPN *a, NM_f_params_t param) which must be define
/// along with other args (params) before including this file
NM_info_t nealder_mead(FPN *a, NM_f_params_t param) {

  FPN V[K][K]; // vertices
  FPN vc[K];   // vertex costs
  FPN vcs[K];  // vc sorted

  FPN c[K];
  FPN r[K];
  FPN t[K];

  int stack[K];
  int ids[K];

  int wi, swi;              // worst and second worst indices
  FPN rc, ec, cco, cci, tf; // fevals and temporary FPN (for casting)

  // Start on the edges of the full simplex
  for (int i = 0; i < K; i++) {
    for (int j = 0; j < K; j++) {
      if (i == j) {
        V[i][j] = FONE;
      } else {
        V[i][j] = FZERO;
      }
    }
  }

  FPN mins[K];
  FPN maxs[K];

  int i = -1;
  // We know c and r will be rewritten when used so we could temporarily use
  // them here (as mins/maxs)
  while (i < N_NM_iter && melbb(V, mins, maxs) > NM_TOL) {
    i++;

    // vertex costs
    for (int k = 0; k < K; k++) {
      for (int j = 0; j < K; j++) {
        t[j] = V[k][j];
      }
      vc[k] = NM_func(t, param);
    }

    quickSortIterative(vc, ids, stack, vcs, K);

    wi = ids[K - 1];
    swi = ids[K - 2];

    // centroid of all but worst
    tf = K - 1;
    for (int j = 0; j < K; j++) {
      c[j] = 0;
      for (int k = 0; k < K - 1; k++) {
        c[j] += V[ids[k]][j] / tf;
      }
    }

    // reflected point
    conv_add(c, V[wi], r, -NM_ALPHA);
    rc = NM_func(r, param);

    if ((vc[ids[0]] <= rc) &&
        (rc < vc[swi])) // cost of reflected in between best and worst vertices,
                        // reflect
    {
      copy_vec(V[wi], r);
      continue;
    }

    if ((vc[ids[0]] >
         rc)) // cost of reflected is below (better than) best vertex,
    {         // expand (try to go further in this direction)
      conv_add(r, c, t,
               -NM_GAMMA); // GAMMA is not used in the usual sense in this
                           // implementation (usual GAMMA = 1 + our GAMMA)
      ec = NM_func(t, param);
      if (ec < rc) {
        copy_vec(V[wi], t);
      } else {
        copy_vec(V[wi], r);
      }
      continue;
    }

    // if we reach here, rc > vc[swi] (worse than second worse)
    if (rc < vc[wi]) // still better than worst
    {
      conv_add(c, r, t, NM_RHO); // "contract outside" (try a point closer to
                                 // center wrt reflected)
      cco = NM_func(t, param);
      if (cco < rc) {
        copy_vec(V[wi], t);
        continue;
      } // fall through to shrink step
    } else {
      conv_add(
          c, V[wi], t,
          NM_RHO); // "contract inside" (try a point closer to center wrt worst)
      cci = NM_func(t, param);
      if (cci < vc[wi]) {
        copy_vec(V[wi], t);
        continue;
      } // fall through to shrink step
    }

    // When all else fails, shrink
    for (int k = 1; k < K; k++) {
      conv_add(V[ids[0]], V[ids[k]], t, NM_SIGMA); // Pull vertex towards best
      copy_vec(V[ids[k]], t);
      // V[ids[k]] could also be manipulated in place
    }
  }

  tf = K;
  // return centroid of final simplex
  for (int j = 0; j < K; j++) {
    a[j] = FZERO;
    for (int k = 0; k < K; k++) {
      a[j] += V[k][j] / tf;
    }
  }

  return (NM_info_t){i, melbb(V, mins, maxs)};
}
