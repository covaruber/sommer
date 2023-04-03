test_that("Interactions with * work properly", {
    data(DT_yatesoats)
    DT <- DT_yatesoats

    m1 <- mmer(fixed=Y ~ V+N+V:N, random = ~ B + B:MP, rcov=~units, data = DT, verbose = F)
    m2 <- mmer(fixed=Y ~ V*N, random = ~ B + B:MP, rcov=~units, data = DT, verbose = F)

    Dt <- m2$Dtable; Dt
    # first fixed effect include and average
    Dt[1,"include"] = TRUE
    Dt[1,"average"] = TRUE
    # second fixed effect include and average
    Dt[2,"include"] = TRUE
    Dt[2,"average"] = TRUE
    # third fixed effect include and average
    Dt[3,"include"] = TRUE
    Dt[3,"average"] = TRUE
    Dt[4,"include"] = TRUE
    Dt[4,"average"] = TRUE
    Dt

    p3 <- predict.mmer(object=m1, D = "V:N", verbose = F, Dtable = Dt)
    p4 <- predict.mmer(object=m2, D = "V:N", verbose = F, Dtable = Dt)
    # expect_equal(p1$pvals$predicted.value, p2$pvals$predicted.value)
    # expect_equal(p1$hypertable$average, p2$hypertable$average)
    expect_equal(p3$pvals$predicted.value, p4$pvals$predicted.value)
    expect_equal(p3$hypertable$average, p4$hypertable$average)
})
