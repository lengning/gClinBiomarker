context("CoxTab")

test_that("Hazard Ratio", {
    data(input)
    res <- CoxTab(data=input,
                  tte="OS",
                  cens="OS.event",
                  var=c('Sex',"Country","Age"))

    expect_that(as.numeric(res["Sex (M/F)", "CI.low"]) <= as.numeric(res["Sex (M/F)", "CI.high"]), is_true())
    expect_that(as.numeric(res["Sex (M/F)", "HR"]) <= as.numeric(res["Sex (M/F)", "CI.high"]), is_true())
    expect_that(as.numeric(res["Sex (M/F)", "HR"]) >= as.numeric(res["Sex (M/F)", "CI.low"]), is_true())
})
