context('utils-io: sink_to_temp')

test_that('sink_to_temp hides output of cat statements', {
  expect_silent(sink_to_temp(cat('test')))
})

test_that('sink_to_temp still outputs warning messages', {
  expect_warning(sink_to_temp(warning('test')), 'test')
})



context('utils-io: print_to_string')

test_that('print_to_string coerces many output lines to a single line', {
  expect_equal(
    length(print_to_string({ print('a'); print('b'); 'c' })),
    1
  )

  expect_equal(
    length(print_to_string(Sys.info())),
    1
  )
})
