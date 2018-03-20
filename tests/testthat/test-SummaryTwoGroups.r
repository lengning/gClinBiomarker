context("SummaryTwoGroups")

test_that("Confidence Interval", {
    data(input)
    res <- SummaryTwoGroups(outcome.var = input$OS,
                       treatment.var = input$Arm,
                       placebo.code = "CTRL",
                       active.code = "TRT",
                       outcome.class = "continuous")

    expect_that(res["Lower"] <= res["Upper"], is_true())
    expect_that(res["P"] <= 1 & res["P"] >= 0, is_true())
    expect_that(res["Effect.Size"] >= res["Lower"] & res["Effect.Size"] <= res["Upper"], is_true())
})
