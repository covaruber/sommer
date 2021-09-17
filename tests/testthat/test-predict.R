test_that("Interactions with * work properly", {
    data(DT_yatesoats)
    DT <- DT_yatesoats

    m1 <- mmer(fixed=Y ~ V+N+V:N, random = ~ B + B:MP, rcov=~units, data = DT, verbose = F)
    m2 <- mmer(fixed=Y ~ V*N, random = ~ B + B:MP, rcov=~units, data = DT, verbose = F)

    p1 <- predict.mmer(m1, classify = "N", verbose = F)
    p2 <- predict.mmer(m2, classify = "N", verbose = F)
    p3 <- predict.mmer(m1, classify = "V:N", verbose = F)
    p4 <- predict.mmer(m2, classify = "V:N", verbose = F)
    expect_equal(p1$pvals$predicted.value, p2$pvals$predicted.value)
    expect_equal(p1$hypertable$average, p2$hypertable$average)
    expect_equal(p3$pvals$predicted.value, p4$pvals$predicted.value)
    expect_equal(p3$hypertable$average, p4$hypertable$average)
})
