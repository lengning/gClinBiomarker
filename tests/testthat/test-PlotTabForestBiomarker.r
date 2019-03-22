context("PlotTabForestBiomarker")

test_that("Confidence Interval", {
    data(input)
    res <- PlotTabForestBiomarker(data=input,
                                  outcome.class="survival",
                                  outcome.var=c("PFS","PFS.event"),
                                  trt="Arm",
                                  var="KRAS.exprs",
                                  var.class="numeric",
                                  var.name="KRAS exprs",
                                  percentile.cutoff=c(.25,.5,.75),
                                  numerical.cutoff=NULL,
                                  greater=TRUE,
                                  less=TRUE,
                                  within.bin=FALSE,
                                  show.itt=TRUE,
                                  show.bep=TRUE,
                                  only.stat=TRUE)

    # All TRT
    LL <- as.numeric(substr(res[3, 6], 1, 4))
    UL <- as.numeric(substr(res[3, 6], 8, 11))
    estimate <- as.numeric(substr(res[3, 6], 8, 11))
    expect_that(LL <= estimate, is_true())
    expect_that(UL >= estimate, is_true())

    # BEP TRT
    LL <- as.numeric(substr(res[5, 6], 1, 4))
    UL <- as.numeric(substr(res[5, 6], 8, 11))
    estimate <- as.numeric(substr(res[5, 6], 8, 11))
    expect_that(LL <= estimate, is_true())
    expect_that(UL >= estimate, is_true())
})
