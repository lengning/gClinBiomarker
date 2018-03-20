context('utils-data-handling: augment_predict')

# test_that('augment_predict: test that deprecation warning suppression is still necessary', {
#   # There's some code in augment_predict to suppress current deprecation
#   # warnings in the broom package. This test will fail when those deprecations
#   # are resolved in broom. When that happens, the deprecation warning
#   # suppression code should be removed from augment_predict.
#   expect_warning(
#     broom::augment(lm(mpg ~ carb + am + hp + wt, mtcars), mtcars),
#     'Deprecated.*purrr::possibly()'
#   )
# })

test_that('augment_predict produces dataframe of model fit terms', {
  expect_equal(
    augment_predict(mtcars, model = lm, model.formula = mpg ~ carb + am + hp + wt),
    suppressWarnings(broom::augment(lm(mpg ~ carb + am + hp + wt, mtcars), mtcars))
  )
})

test_that('augment_predict produces grouped data frame when multiple models requested', {
  expect_true(
    'grouped_df' %in%
    class(augment_predict(mtcars, model = lm, model.per = ~ am,
      model.formula = mpg ~ carb + am + hp + wt))
  )

  expect_equal(
    attr(augment_predict(mtcars, model = lm, model.per = ~ am,
      model.formula = mpg ~ carb + am + hp + wt), 'vars'),
    'am'
  )
})



context('utils-data-handling: copy_variable_attributes')

test_that('ensure that variable attributes migrate properly between dataframes', {
  test_mtcars <- mtcars
  attr(test_mtcars[['carb']], 'attr1') <- 'test1'
  attr(test_mtcars[['carb']], 'attr2') <- 'test2'
  attr(test_mtcars[['cyl']], 'attr1') <- 'test1'

  expect_equal(
    copy_variable_attributes(mtcars, test_mtcars),
    test_mtcars
  )
})



context('utils-data-handling: summarize_unary_vars')

test_that('summarize_unary_vars creates new group_by dataframe of unary value columns', {
  test_mtcars <- mtcars
  test_mtcars[['carb2']] <- test_mtcars[['carb']] * 2
  test_mtcars[['carb4']] <- test_mtcars[['carb']] * 4

  expect_equal(
    names(summarize_unary_vars(test_mtcars, 'carb')),
    c('carb', 'carb2', 'carb4')
  )

  expect_equal(
    nrow(summarize_unary_vars(test_mtcars, 'carb')),
    length(unique(test_mtcars[['carb']]))
  )
})
