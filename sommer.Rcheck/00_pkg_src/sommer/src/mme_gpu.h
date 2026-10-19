#ifndef SOMMER_MME_GPU_H
#define SOMMER_MME_GPU_H

#include <RcppArmadillo.h>
#include <climits>
#include <limits>
#include <memory>
#include <string>
#include "mme_cuda_runtime.h"

namespace sommer {
namespace gpu {

class DenseBackend {
public:
  explicit DenseBackend(bool enabled, int minimumSize = 256)
    : enabled_(enabled), minimumSize_(minimumSize) {
    if(minimumSize < 1){ Rcpp::stop("gpuMinSize must be positive."); }
#ifdef SOMMER_HAVE_CUDA
    if(enabled_){
      runtime_.reset(new Runtime());
      check(runtime_->cudaFree(nullptr), "initialize CUDA");
      if(runtime_->cublasCreate(&blas_) != CUBLAS_STATUS_SUCCESS){
        Rcpp::stop("Unable to create cuBLAS handle.");
      }
      if(runtime_->cusolverDnCreate(&solver_) != CUSOLVER_STATUS_SUCCESS){
        runtime_->cublasDestroy(blas_);
        blas_ = nullptr;
        Rcpp::stop("Unable to create cuSOLVER handle.");
      }
    }
#else
    if(enabled_){
      Rcpp::stop("GPU kernels require a CUDA-enabled sommer build. Install the CUDA toolkit and reinstall with SOMMER_CUDA_HOME set; use gpu=FALSE for CPU fitting.");
    }
#endif
  }

  ~DenseBackend(){
#ifdef SOMMER_HAVE_CUDA
    if(solver_){ runtime_->cusolverDnDestroy(solver_); }
    if(blas_){ runtime_->cublasDestroy(blas_); }
#endif
  }
  DenseBackend(const DenseBackend &) = delete;
  DenseBackend & operator=(const DenseBackend &) = delete;

  bool cholesky(arma::mat & output, const arma::mat & input,
                const char * layout = "upper") {
    if(!eligible(input.n_rows)){ return arma::chol(output, input, layout); }
#ifdef SOMMER_HAVE_CUDA
    if(input.n_rows != input.n_cols){ Rcpp::stop("GPU Cholesky requires a square matrix."); }
    const int dimension = checkedDimension(input.n_rows);
    const cublasFillMode_t triangle = std::string(layout) == "lower"
      ? CUBLAS_FILL_MODE_LOWER : CUBLAS_FILL_MODE_UPPER;
    Buffer<double> matrix(*runtime_, input.n_elem);
    upload(matrix.get(), input.memptr(), input.n_elem);
    int workspaceSize = 0;
    checkSolver(runtime_->cusolverDnDpotrf_bufferSize(solver_, triangle, dimension,
      matrix.get(), dimension, &workspaceSize), "Cholesky workspace query");
    Buffer<double> workspace(*runtime_, static_cast<std::size_t>(workspaceSize));
    Buffer<int> info(*runtime_, 1);
    checkSolver(runtime_->cusolverDnDpotrf(solver_, triangle, dimension, matrix.get(),
      dimension, workspace.get(), workspaceSize, info.get()), "Cholesky factorization");
    int status = 0;
    download(&status, info.get(), 1);
    if(status < 0){ Rcpp::stop("cuSOLVER Cholesky received an invalid argument."); }
    if(status > 0){ output.reset(); return false; }
    output.set_size(input.n_rows, input.n_cols);
    download(output.memptr(), matrix.get(), input.n_elem);
    output = triangle == CUBLAS_FILL_MODE_LOWER
      ? arma::mat(arma::trimatl(output)) : arma::mat(arma::trimatu(output));
    ++choleskyCalls;
    return true;
#else
    return false;
#endif
  }

  template<typename MatrixExpression, typename RhsExpression>
  bool triangularSolve(arma::mat & output,
      const arma::Op<MatrixExpression, arma::op_trimat> & triangle,
      const RhsExpression & rhs, const arma::solve_opts::opts & options) {
    if(!enabled_){ return arma::solve(output, triangle, rhs, options); }
    const arma::mat matrix(triangle);
    const arma::mat right(rhs);
    if(!eligible(matrix.n_rows) || right.n_cols == 0){
      return arma::solve(output, triangle, rhs, options);
    }
#ifdef SOMMER_HAVE_CUDA
    if(matrix.n_rows != matrix.n_cols || right.n_rows != matrix.n_rows){
      Rcpp::stop("GPU triangular solve dimensions do not match.");
    }
    if(!matrix.is_finite() || !right.is_finite() || arma::any(matrix.diag() == 0.0)){
      return false;
    }
    const int dimension = checkedDimension(matrix.n_rows);
    const int columns = checkedDimension(right.n_cols);
    Buffer<double> deviceMatrix(*runtime_, matrix.n_elem);
    Buffer<double> deviceRight(*runtime_, right.n_elem);
    upload(deviceMatrix.get(), matrix.memptr(), matrix.n_elem);
    upload(deviceRight.get(), right.memptr(), right.n_elem);
    const double alpha = 1.0;
    const cublasFillMode_t fill = triangle.aux_uword_a == 0
      ? CUBLAS_FILL_MODE_UPPER : CUBLAS_FILL_MODE_LOWER;
    if(runtime_->cublasDtrsm(blas_, CUBLAS_SIDE_LEFT, fill, CUBLAS_OP_N,
        CUBLAS_DIAG_NON_UNIT, dimension, columns, &alpha, deviceMatrix.get(),
        dimension, deviceRight.get(), dimension) != CUBLAS_STATUS_SUCCESS){
      Rcpp::stop("cuBLAS triangular solve failed.");
    }
    output.set_size(right.n_rows, right.n_cols);
    download(output.memptr(), deviceRight.get(), right.n_elem);
    ++triangularSolveCalls;
    return output.is_finite();
#else
    return false;
#endif
  }

  std::size_t choleskyCalls = 0;
  std::size_t triangularSolveCalls = 0;
  std::size_t productCalls = 0;
  bool enabled() const { return enabled_; }

  arma::mat product(const arma::mat & left, const arma::mat & right,
                    bool transposeLeft = false, bool transposeRight = false){
    const arma::uword rows = transposeLeft ? left.n_cols : left.n_rows;
    const arma::uword inner = transposeLeft ? left.n_rows : left.n_cols;
    const arma::uword columns = transposeRight ? right.n_rows : right.n_cols;
    const arma::uword rightInner = transposeRight ? right.n_cols : right.n_rows;
    if(inner != rightInner){ Rcpp::stop("GPU matrix product dimensions do not match."); }
    if(!eligible(rows) || !eligible(columns) || !eligible(inner)){
      if(transposeLeft && transposeRight){ return left.t() * right.t(); }
      if(transposeLeft){ return left.t() * right; }
      if(transposeRight){ return left * right.t(); }
      return left * right;
    }
#ifdef SOMMER_HAVE_CUDA
    arma::mat output(rows, columns);
    Buffer<double> deviceLeft(*runtime_, left.n_elem);
    Buffer<double> deviceRight(*runtime_, right.n_elem);
    Buffer<double> deviceOutput(*runtime_, output.n_elem);
    upload(deviceLeft.get(), left.memptr(), left.n_elem);
    upload(deviceRight.get(), right.memptr(), right.n_elem);
    const double alpha = 1.0;
    const double beta = 0.0;
    if(runtime_->cublasDgemm(blas_, transposeLeft ? CUBLAS_OP_T : CUBLAS_OP_N,
        transposeRight ? CUBLAS_OP_T : CUBLAS_OP_N, checkedDimension(rows),
        checkedDimension(columns), checkedDimension(inner), &alpha,
        deviceLeft.get(), checkedDimension(left.n_rows), deviceRight.get(),
        checkedDimension(right.n_rows), &beta, deviceOutput.get(), checkedDimension(rows))
        != CUBLAS_STATUS_SUCCESS){
      Rcpp::stop("cuBLAS matrix product failed.");
    }
    download(output.memptr(), deviceOutput.get(), output.n_elem);
    ++productCalls;
    return output;
#else
    return arma::mat();
#endif
  }

private:
  bool enabled_;
  int minimumSize_;
  bool eligible(arma::uword size) const {
    return enabled_ && size >= static_cast<arma::uword>(minimumSize_);
  }
#ifdef SOMMER_HAVE_CUDA
  std::unique_ptr<Runtime> runtime_;
  cublasHandle_t blas_ = nullptr;
  cusolverDnHandle_t solver_ = nullptr;

  static int checkedDimension(arma::uword size){
    if(size > INT_MAX){ Rcpp::stop("GPU matrix dimension exceeds CUDA integer limits."); }
    return static_cast<int>(size);
  }
  void check(cudaError_t status, const char * operation){
    if(status != cudaSuccess){
      Rcpp::stop(std::string("CUDA ") + operation + ": " + runtime_->cudaGetErrorString(status));
    }
  }
  static void checkSolver(cusolverStatus_t status, const char * operation){
    if(status != CUSOLVER_STATUS_SUCCESS){
      Rcpp::stop(std::string("cuSOLVER failed during ") + operation + ".");
    }
  }
  template<typename Value> class Buffer {
  public:
    explicit Buffer(Runtime & runtime, std::size_t count) : runtime_(runtime){
      if(count > std::numeric_limits<std::size_t>::max() / sizeof(Value)){
        Rcpp::stop("GPU allocation size overflow.");
      }
      const std::size_t bytes = count * sizeof(Value);
      const cudaError_t status = runtime_.cudaMalloc(reinterpret_cast<void **>(&pointer_), bytes);
      if(status != cudaSuccess){
        Rcpp::stop("GPU allocation of %.0f bytes failed: %s. Use a GPU with more free VRAM or reduce model size; no CPU fallback was attempted.",
          static_cast<double>(bytes), runtime_.cudaGetErrorString(status));
      }
    }
    ~Buffer(){ if(pointer_){ runtime_.cudaFree(pointer_); } }
    Buffer(const Buffer &) = delete;
    Buffer & operator=(const Buffer &) = delete;
    Value * get() const { return pointer_; }
  private:
    Runtime & runtime_;
    Value * pointer_ = nullptr;
  };
  template<typename Value>
  void upload(Value * device, const Value * host, std::size_t count){
    check(runtime_->cudaMemcpy(device, host, count * sizeof(Value), cudaMemcpyHostToDevice), "upload");
  }
  template<typename Value>
  void download(Value * host, const Value * device, std::size_t count){
    check(runtime_->cudaMemcpy(host, device, count * sizeof(Value), cudaMemcpyDeviceToHost), "download");
  }
#endif
};

}
}
#endif
