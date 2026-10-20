test_that("row-wise spatial products preserve ordering and singleton dimensions", {
  for(dimensions in list(c(1L, 2L, 3L), c(2L, 1L, 1L), c(4L, 2L, 3L))){
    n <- dimensions[1L]
    A <- matrix(seq_len(n * dimensions[2L]), n)
    B <- matrix(seq_len(n * dimensions[3L]), n)
    expected <- matrix(unlist(lapply(seq_len(n),
      function(i) kronecker(A[i,], B[i,]))), nrow=n, byrow=TRUE)
    expect_equal(sommer:::.spatial_row_tensor(A, B), expected,
                 ignore_attr=TRUE)
    expect_equal(dim(sommer:::.spatial_row_tensor(A, B)), dim(expected))
  }
})

test_that("spatial spline designs preserve repeated observations", {
  data <- expand.grid(col=seq_len(6L), row=seq_len(5L))
  data <- data[-c(2L, 8L, 19L),]
  repeated <- rbind(data, data[1L,])
  spline <- spl2Dc(repeated$col, repeated$row, nsegments=c(6L, 4L))
  design <- as.matrix(spline$Z[[1L]])
  expect_equal(dim(design), c(nrow(repeated), nrow(spline$Gu)))
  expect_equal(design[1L,], design[nrow(design),])
  expect_true(all(is.finite(design)))
  for(method in c("Lee", "Wood")){
    basis <- spl2Dmats("col", "row", data=data,
                      nsegments=c(6L, 4L), method=method)
    expect_equal(nrow(basis$all), nrow(data))
    components <- do.call(cbind, basis[seq_len(5L)])
    expect_equal(basis$all, components[, colnames(basis$all)], ignore_attr=TRUE)
    expect_true(all(is.finite(basis$all)))
  }
})
