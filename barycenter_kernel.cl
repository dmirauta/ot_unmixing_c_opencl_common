#ifndef EXTERNAL_CONCAT
 #include "shared_defines.h"
 #include "dists.c"
 #include "lb.c"
#endif

__kernel void multibary_pass(__global FPN *XC_g,
                             __global FPN *WC_g,
                             __global FPN *λv_g,
                                      int  skip)
{
    int i = get_global_id(0);
        
    FPN xb_local[N];
    
    int il = 2*i*skip;
    int iu = il+skip;
    int il_start = il*N; //start of lower item in pair
    int iu_start = iu*N;
    
    FPN λ = λv_g[il]/(λv_g[il] + λv_g[iu]+EPSILON);

    int path_a[N * 2 - 1];
    int path_b[N * 2 - 1];
    FPN transfer[N * 2 - 1];
    TraceData_t td = {transfer, path_a, path_b};
    Signature_t sa = {
        &XC_g[il_start],
        &WC_g[il_start],
        N
    };
    Signature_t sb = {
        &XC_g[iu_start],
        &WC_g[iu_start],
        N
    };

    left_bary(sa,
              sb,
              xb_local,
              λ,
              td,
              EPSILON);
    
    // write bary into lower signature of each pair
    for(int j=0; j<N; j++)
    {
        XC_g[il_start + j] = xb_local[j];
    }
    λv_g[il] += λv_g[iu];
    //λv_g[iu]=0; // not necessary, but marks no longer used positions
    
}
