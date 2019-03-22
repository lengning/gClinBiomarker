context("SummaryVars")

test_that("Size of data", {
    data(input)
    res <- SummaryVars(data=input,
                       trt='Arm',
                       subgroup='BEP',
                       var=c('Age','Sex'),
                       var.class=c("numeric","categorical"))

    expect_equal(as.numeric(res["N", "All(CTRL)"]) + as.numeric(res["N", "All(TRT)"]), nrow(input))
    expect_equal(as.numeric(res["Total (non-NA)", "All(CTRL)"]) + as.numeric(res["Total (non-NA)", "All(TRT)"]), nrow(input))

    expect_equal(as.numeric(res["Total (non-NA)", "All(CTRL)"]), as.numeric(res["N", "All(CTRL)"]))
    expect_equal(as.numeric(res["Total (non-NA)", "BEP(CTRL)"]), as.numeric(res["N", "BEP(CTRL)"]))
    expect_equal(as.numeric(res["Total (non-NA)", "All(TRT)"]), as.numeric(res["N", "All(TRT)"]))
    expect_equal(as.numeric(res["Total (non-NA)", "BEP(TRT)"]), as.numeric(res["N", "BEP(TRT)"]))
})
