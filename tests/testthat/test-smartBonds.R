test_that("smartBonds() returns NULL if smiles input is NULL or NA", {
  expect_type(smartBonds(NA, "C=C"), "NULL")
  expect_type(smartBonds(NULL, "C=C"), "NULL")
})
