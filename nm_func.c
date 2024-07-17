
typedef struct NM_f_params {
  const FPN *xc;
  const FPN *wc;
  GLB const FPN *XBfs; // passed through flat and sorted
  GLB const FPN *WBfs; // passed through flat and sorted
  GLB const int *karr; // basis correspondence
  FPN *cts;
} NM_f_params_t;

FPN NM_func(FPN *a, NM_f_params_t param) {
  return segment_cost_aug(param.xc, param.wc, param.XBfs, param.WBfs, a,
                          param.cts, param.karr);
}
