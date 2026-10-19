#ifndef SOMMER_MME_CUDA_RUNTIME_H
#define SOMMER_MME_CUDA_RUNTIME_H

#include <cstdlib>
#include <stdexcept>
#include <string>
#include <vector>
#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <dlfcn.h>
#endif
#ifdef SOMMER_HAVE_CUDA
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#include <cusolverDn.h>
#endif

namespace sommer {
namespace gpu {

class Library {
public:
  explicit Library(const std::vector<std::string> & names){
    std::vector<std::string> directories{""};
    const char * home = std::getenv("SOMMER_CUDA_HOME");
    if(home){
      directories.push_back(std::string(home) + "/lib64/");
      directories.push_back(std::string(home) + "/lib/");
      directories.push_back(std::string(home) + "/bin/");
    }
    directories.push_back("/usr/local/cuda/lib64/");
    for(const std::string & directory : directories){
      for(const std::string & name : names){
#ifdef _WIN32
        handle_ = LoadLibraryA((directory + name).c_str());
#else
        handle_ = dlopen((directory + name).c_str(), RTLD_NOW | RTLD_LOCAL);
#endif
        if(handle_){ return; }
      }
    }
    throw std::runtime_error("Cannot load CUDA library " + names.front() +
      ". Install CUDA runtime, cuBLAS and cuSOLVER and expose their library directory through SOMMER_CUDA_HOME or the system loader path.");
  }
  ~Library(){
#ifdef _WIN32
    if(handle_){ FreeLibrary(handle_); }
#else
    if(handle_){ dlclose(handle_); }
#endif
  }
  Library(const Library &) = delete;
  Library & operator=(const Library &) = delete;
  template<typename Function> Function symbol(const char * name){
#ifdef _WIN32
    const auto address = GetProcAddress(handle_, name);
#else
    void * address = dlsym(handle_, name);
#endif
    if(!address){ throw std::runtime_error(std::string("Missing CUDA symbol ") + name + ". Install compatible CUDA libraries."); }
    return reinterpret_cast<Function>(address);
  }
private:
#ifdef _WIN32
  HMODULE handle_ = nullptr;
#else
  void * handle_ = nullptr;
#endif
};

struct Runtime {
#ifdef _WIN32
  Library runtime{{"cudart64_12.dll", "cudart64_13.dll", "cudart64_110.dll"}};
  Library blas{{"cublas64_12.dll", "cublas64_13.dll", "cublas64_11.dll"}};
  Library solver{{"cusolver64_11.dll", "cusolver64_12.dll", "cusolver64_10.dll"}};
#else
  Library runtime{{"libcudart.so", "libcudart.so.12", "libcudart.so.13", "libcudart.so.11.0"}};
  Library blas{{"libcublas.so", "libcublas.so.12", "libcublas.so.13", "libcublas.so.11"}};
  Library solver{{"libcusolver.so", "libcusolver.so.11", "libcusolver.so.12", "libcusolver.so.10"}};
#endif
#ifdef SOMMER_HAVE_CUDA
#define SOMMER_GPU_STRING_INNER(value) #value
#define SOMMER_GPU_STRING(value) SOMMER_GPU_STRING_INNER(value)
#define SOMMER_GPU_SYMBOL(library, name) decltype(&::name) name = library.symbol<decltype(&::name)>(SOMMER_GPU_STRING(name))
  SOMMER_GPU_SYMBOL(runtime, cudaFree);
  SOMMER_GPU_SYMBOL(runtime, cudaMalloc);
  SOMMER_GPU_SYMBOL(runtime, cudaMemcpy);
  SOMMER_GPU_SYMBOL(runtime, cudaGetErrorString);
  SOMMER_GPU_SYMBOL(runtime, cudaGetDeviceCount);
  SOMMER_GPU_SYMBOL(blas, cublasCreate);
  SOMMER_GPU_SYMBOL(blas, cublasDestroy);
  SOMMER_GPU_SYMBOL(blas, cublasDtrsm);
  SOMMER_GPU_SYMBOL(blas, cublasDgemm);
  SOMMER_GPU_SYMBOL(solver, cusolverDnCreate);
  SOMMER_GPU_SYMBOL(solver, cusolverDnDestroy);
  SOMMER_GPU_SYMBOL(solver, cusolverDnDpotrf_bufferSize);
  SOMMER_GPU_SYMBOL(solver, cusolverDnDpotrf);
#undef SOMMER_GPU_SYMBOL
#undef SOMMER_GPU_STRING
#undef SOMMER_GPU_STRING_INNER
#endif
};

}
}
#endif
