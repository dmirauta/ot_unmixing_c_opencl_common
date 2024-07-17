#ifndef EXTERNAL_CONCAT
#include "shared_defines.h"
#include "segment_cost.c"
#include "nm_func.c"
#include "nealder_mead.c"
#endif

__kernel void multisegment( __global const FPN *XC_g,
                            __global const FPN *WC_g,
                            __global const FPN *XBfs_g,
                            __global const FPN *WBfs_g,
                            __global const int *karr_g,
                            __global       FPN *Bta_g,
                            __global       FPN *Wa_g,
                            __global       int *Ia_g 
                            )
{
    int i = get_global_id(0);
        
    if(i<M)
    {
    FPN xc_local[N];
    FPN wc_local[N];
    FPN a_local[K];
    FPN Bt[K*N];
    
    // Put into local memory
    for(int j=0; j<N; j++)
    {
        xc_local[j] = XC_g[i*N + j];
        wc_local[j] = WC_g[i*N + j];
    }
    
    FPN cts[K*N];
    NM_f_params_t param = {xc_local,
                           wc_local,
                           XBfs_g,
                           WBfs_g,
                           karr_g,
                           cts};

    NM_info_t nm_info = nealder_mead(a_local, param);
    Ia_g[i] = nm_info.iter; // debug
    
    FPN a_sum = EPSILON;
    for(int k=0; k<K; k++)
    {
        a_sum += a_local[k]; 
    }
    for(int k=0; k<K; k++)
    {
        a_local[k] /= a_sum;
    }
    
    Wa_g[i] = segment_cost_rc(xc_local,
                                wc_local,
                                XBfs_g,
                                WBfs_g, // passed through flat and sorted
                                a_local,
                                cts,
                                karr_g,
                                Bt);
    //Wa_g[i] = nm_info.diam; // debug
    
    // Export to global
    for(int j=0; j<K; j++)
    {
        
        for(int k=0; k<N; k++)
        {
            Bta_g[i*N*K + j*N + k] = Bt[j*N + k];
        }
    }
    
    }
    
}
