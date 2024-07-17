#ifndef UNMIX_SHARED_DEFINES

#ifdef USE_FLOAT
typedef float FPN;
#define FZERO 0.0f
#define FONE 1.0f
#else
typedef double FPN;
#define FZERO 0.0
#define FONE 1.0
#endif

#ifdef __OPENCL_VERSION__
#define GLB __global
#else
#define GLB // ignore CL memory scopes when compiled on host
#endif

// NOTE: a data wrapper for more concise function signatures
typedef struct Signature {
  /// positions
  GLB FPN const *x;
  /// weights
  GLB FPN const *w;
  /// size
  int n;
} Signature_t;

/// defines a monotone transport plan.
/// supplied externally even when just required locally to leave assigment
/// details up to caller
typedef struct TraceData {
  /// amount transferred between corresponding positions
  FPN *transfer;
  /// indeces into positions of marginal a
  int *path_a;
  /// indeces into positions of marginal b
  int *path_b;
} TraceData_t;

#define UNMIX_SHARED_DEFINES
#endif
