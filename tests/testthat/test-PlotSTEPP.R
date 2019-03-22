context("PlotSTEPP")

test_that("Hazard Ratio", {
    data(input)
    res <- PlotSTEPP(data = input,
                     outcome.var = c("PFS", "PFS.event"),
                     outcome.class = "survival",
                     trt = "Arm",
                     var = "KRAS.exprs",
                     covariate = "Sex",
                     strata = "Age",
                     placebo.code = "CTRL",
                     active.code = "TRT")

    expect_that(as.numeric(res[1, "Hazard Ratio"]) >= as.numeric(res[1, "CI Lower"]), is_true())
    expect_that(as.numeric(res[1, "Hazard Ratio"]) <= as.numeric(res[1, "CI Upper"]), is_true())
    expect_that(as.numeric(res[1, "CI Lower"]) <= as.numeric(res[1, "CI Upper"]), is_true())

    for (i in 1:16) {
        expect_that(as.numeric(res[i, "Window Left"]) < as.numeric(res[i, "Window Right"]), is_true())
    }
})
