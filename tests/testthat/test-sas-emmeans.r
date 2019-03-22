## setup
cols <- c('carb', 'hp', 'wt')

test_mtcars <- mtcars
test_mtcars$carb <- as.factor(test_mtcars$carb)
test_mtcars$cyl <- as.factor(test_mtcars$cyl)



context("sas.emmeans() helper functions: numeric_means")

test_that('numeric_means returns a list', {
  expect_equal(class(numeric_means(mtcars)), 'list')
})

test_that('numeric_means defaults to calculating means of all numeric vars', {
  expect_equal(
    names(numeric_means(mtcars)),
    names(mtcars)
  )
  expect_equal(
    names(numeric_means(test_mtcars)),
    setdiff(names(mtcars), c('carb', 'cyl'))
  )
})

test_that('numeric_means can be set to calculate only a subset of vars', {
  expect_equal(
    names(numeric_means(mtcars, cols)),
    cols
  )
  expect_equal(
    names(numeric_means(test_mtcars, cols)),
    setdiff(cols, 'carb')
  )
})

test_that('numeric_means verbose info is printed when variables omitted', {
  expect_message(
    gClinBiomarker:::numeric_means(test_mtcars, verbose = TRUE),
    'carb'
  )
})



context("sas.emmeans() helper functions: append_weights")

test_that('append_weights returns weights of "proportional" when no levels produced', {
  emmeans_args <- list(
    lm(mpg ~ carb + wt + hp, data = mtcars),
    data = mtcars,
    specs = ~ carb
  )

  expect_equal(
    append_weights(emmeans_args, data = mtcars)$weights,
    'proportional'
  )
})

test_that('append_weights calculates weights over factor variables in specs', {
  emmeans_args <- list(
    lm(mpg ~ carb + wt + hp + vs + am, data = test_mtcars),
    data = test_mtcars,
    specs = pairwise ~ vs + am
  )

  expect_equal(
    append_weights(emmeans_args)$weights,
    as.data.frame(table(mtcars$carb))$Freq
  )
})

test_that('append_weights shows output for verbose output', {
  emmeans_args <- list(
    lm(mpg ~ carb + wt + hp + vs + am, data = test_mtcars),
    data = test_mtcars,
    specs = pairwise ~ vs + am
  )

  expect_message(
    append_weights(emmeans_args, verbose = TRUE),
    ': carb'
  )

  expect_message(
    append_weights(emmeans_args, verbose = TRUE),
    'wgt'
  )
})



context("sas.emmeans() helper functions: coerce_from_factor")

coerce_list <- list(carb = 'character', cyl = 'numeric')

test_that('coerce_from_factor appropriately coerces variables', {
  expect_equal(
    class(coerce_from_factor(test_mtcars, coerce_list)$carb),
    coerce_list$carb
  )

  expect_equal(
    class(coerce_from_factor(test_mtcars, coerce_list)$cyl),
    coerce_list$cyl
  )
})

test_that('coerce_from_factor old vars preserved correctly', {
  expect_equal(
    class(coerce_from_factor(test_mtcars, coerce_list)$carb.old),
    'factor'
  )

  expect_equal(
    class(coerce_from_factor(test_mtcars,
      coerce_list, suffix = '.test')$carb.test),
    'factor'
  )

  expect_equal(
    coerce_from_factor(test_mtcars, coerce_list, preserve = FALSE)$carb.old,
    NULL
  )
})

test_that('coerce_from_factor with no class list retains data', {
  expect_equal(
    coerce_from_factor(mtcars, list()),
    mtcars
  )
})



context("sas.emmeans() helper functions: get_model_formula")

test_func <- function(model) {
  mydata <- mtcars
  model(mpg ~ carb + wt + hp, mydata)
}
# create some models in another environment
model_w_env <- test_func(lm) # lm preserves model initialization environment
model_wo_env <- test_func(nlme::gls) # gls does not

test_that('get_model_formula successfully strips model when it\'s available in model object', {
  expect_false({
    identical(
      environment(),
      attr(terms(as.formula(get_model_formula(model_w_env))), ".Environment")
    )
  })

  expect_true({
    identical(
      environment(),
      attr(terms(as.formula(get_model_formula(model_wo_env))), ".Environment")
    )
  })
})

test_that('get_model_formula fails appropriately when invalid calls are passed', {
  expect_error(
    get_model_formula(list()),
    'Model object does not contain any arguments.'
  )

  expect_error(
    gClinBiomarker:::get_model_formula(list(call = list('function_name', a = 1))),
    'Attempted to coerce first argument'
  )
})



context("sas.emmeans() helper functions: append_split_terms")

test_that('append_split_terms have appropriate list elements', {
  t <- append_split_terms(a ~ b + c + d)
  expect_equal(t$resp, 'a')
  expect_equal(t$covs, c('b', 'c', 'd'))
  expect_equal(t$vars, c('a', 'b', 'c', 'd'))

  t <- append_split_terms(~ b + c + d)
  expect_equal(t$resp, c())
  expect_equal(t$covs, c('b', 'c', 'd'))
  expect_equal(t$vars, c('b', 'c', 'd'))
})



context("sas.emmeans() helper functions: clean_emmeans_specs")

test_that('clean_emmeans_specs appropriately adds response var', {
  f <- ~ a + b
  expect_true(attr(terms(clean_emmeans_specs(f)), 'response') == 1)
  expect_true(terms(clean_emmeans_specs(f))[[2]] == 'pairwise')

  non_formula <- list(pairwise = c('a', 'b'))
  expect_equal(clean_emmeans_specs(non_formula), non_formula)
})



context("sas.emmeans() helper functions: get_emmeans_specs_predictors")

test_that('specs as formula omits response variable', {
  specs <- a ~ b + c
  expect_equal(
    get_emmeans_specs_predictors(specs),
    c('b', 'c')
  )
})

test_that('specs as character returns all', {
  specs <- c('a', 'b', 'c')
  expect_equal(
    get_emmeans_specs_predictors(specs),
    c('a', 'b', 'c')
  )
})

test_that('specs as list notifies of lack of support', {
  specs <- list(pairwise = c('a', 'b'), contrasts = c('c', 'd'))
  expect_error(
    get_emmeans_specs_predictors(specs),
    'Passing list.*unsupported'
  )
})

test_that('specs as any other class throws error', {
  specs <- as.data.frame(list(a = c(1, 2), b = c(3, 4)))
  expect_error(
    get_emmeans_specs_predictors(specs),
    'Unable to parse'
  )
})



context("sas.emmeans()")

test_that('sas.emmeans runs', {
  model <- lm(mpg ~ carb + cyl + wt + am, data = test_mtcars)
  expect_silent(sas.emmeans(model, ~ carb + cyl, quietly = T))
})
