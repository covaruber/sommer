test_that("predict works", {
    data(DT_yatesoats)
    DT <- DT_yatesoats

    m1 <- mmer(fixed=Y ~ V+N+V:N, random = ~ B + B:MP, rcov=~units, data = DT, verbose = F)
    # m2 <- mmer(fixed=Y ~ V*N, random = ~ B + B:MP, rcov=~units, data = DT)

    p1 <- predict.mmer(m1, classify = "N", verbose = F)
    # p2 <- predict.mmer(m2, classify = "N")
    expect_equal(p1$pvals$predicted.value, c(79.38889, 98.88889, 114.22222, 123.38889))
    expect_equal(p1$hypertable$average, c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE))
})
