#include "mme_gpu.h"

// [[Rcpp::export]]
Rcpp::List mme_gpu_runtime_status(){
#ifdef SOMMER_HAVE_CUDA
  const bool compiled = true;
#else
  const bool compiled = false;
#endif
  bool libraries = false;
  bool available = false;
  std::string reason;
  try {
    sommer::gpu::Runtime runtime;
    libraries = true;
#ifdef SOMMER_HAVE_CUDA
    int count = 0;
    const cudaError_t status = runtime.cudaGetDeviceCount(&count);
    if(status != cudaSuccess){
      reason = runtime.cudaGetErrorString(status);
    }else if(count < 1){
      reason = "CUDA reports no accessible NVIDIA devices; check driver installation and CUDA_VISIBLE_DEVICES.";
    }else{
      sommer::gpu::DenseBackend backend(true);
      available = true;
    }
#else
    reason = "CUDA libraries load, but sommer was compiled without CUDA kernels.";
#endif
  } catch(const std::exception & error){
    reason = error.what();
  }
  return Rcpp::List::create(Rcpp::Named("libraries") = libraries,
    Rcpp::Named("compiled") = compiled, Rcpp::Named("available") = available,
    Rcpp::Named("reason") = reason);
}
