#ifndef EXTERNAL_CONCAT
#include "shared_defines.h"
#include "dists.c"
#endif

__kernel void sq_dist_mat(const __global FPN *XC_g,
                          const __global FPN *WC_g,
                                __global FPN *D_g)
{
  int i = get_global_id(0);
  int j = get_global_id(1);
  int M = get_global_size(0);

  int path_a[N * 2 - 1];
  int path_b[N * 2 - 1];
  FPN transfer[N * 2 - 1];
  TraceData_t td = {transfer, path_a, path_b};
  Signature_t sa = {
      &XC_g[i*N],
      &WC_g[i*N],
      N
  };
  Signature_t sb = {
      &XC_g[j*N],
      &WC_g[j*N],
      N
  };
        
  if(i<M) // TODO: redundant check?
  {
    if (j<i) {
        D_g[M*i + j] = W22(sa, sb, td);
    } else {
        D_g[M*i + j] = 0.0f;
    }

  }
  
}

__kernel void sq_dist_mat_from( const __global FPN *XC_g,
                                const __global FPN *WC_g,
                                const __global FPN *XB_g,
                                const __global FPN *WB_g,
                                      __global FPN *D_g)
{
  int i = get_global_id(0);
  int j = get_global_id(1);
  int M = get_global_size(0);
  int K = get_global_size(1);
        
  int path_a[N * 2 - 1];
  int path_b[N * 2 - 1];
  FPN transfer[N * 2 - 1];
  TraceData_t td = {transfer, path_a, path_b};
  Signature_t sa = {
      &XC_g[i*N],
      &WC_g[i*N],
      N
  };
  Signature_t sb = {
      &XB_g[j*N],
      &WB_g[j*N],
      N
  };

  if(i<M && j<K)
  {
    D_g[K*i + j] = W22(sa, sb, td);

  }
  
}
