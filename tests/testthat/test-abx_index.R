testthat::context("Testing abxindx functionality")

testthat::test_that("Testing abxidx on abx_test_df", {

  ##vanco test
  vanco_idx_outcome <- vancomycin_index(abx_test_df)

  testthat::expect_type(vanco_idx_outcome, "double")
  testthat::expect_equal(length(vanco_idx_outcome), 4)
  testthat::expect_equal(vanco_idx_outcome, c(a = log10(0.5/0.5),
                                              b = log10(0.75/0.25),
                                              c = log10(0/1),
                                              d = log10(1/0)))

  ##tetracycline test
  tet_idx_outcome <- tetracycline_index(abx_test_df)

  testthat::expect_type(tet_idx_outcome, "double")
  testthat::expect_equal(length(tet_idx_outcome), 4)
  testthat::expect_equal(tet_idx_outcome, c(a = log10(0.4/0.6),
                                            b = log10(0.25/0.75),
                                            c = log10(0/1),
                                            d = log10(0/1)))

  ##gram positive test
  gram_pos_idx_outcome <- gram_pos_index(abx_test_df)

  testthat::expect_type(gram_pos_idx_outcome, "double")
  testthat::expect_equal(length(gram_pos_idx_outcome), 4)
  testthat::expect_equal(gram_pos_idx_outcome, c(a = log10(0.4/0.6),
                                                 b = log10(0.5/0.5),
                                                 c = log10(0/1),
                                                 d = log10(0/1)))

  ##gram negative test
  gram_neg_idx_outcome <- gram_neg_index(abx_test_df)

  testthat::expect_type(gram_neg_idx_outcome, "double")
  testthat::expect_equal(length(gram_neg_idx_outcome), 4)
  testthat::expect_equal(gram_neg_idx_outcome, c(a = log10(0.6/0.4),
                                                 b = log10(0.5/0.5),
                                                 c = log10(1/0),
                                                 d = log10(1/0)))

  ##anaerobes test
  ana_idx_outcome <- anaerobes_index(abx_test_df)

  testthat::expect_type(ana_idx_outcome, "double")
  testthat::expect_equal(length(ana_idx_outcome), 4)
  testthat::expect_equal(ana_idx_outcome, c(a = log10(0.7/0.3),
                                            b = log10(0.75/0.25),
                                            c = log10(1/0),
                                            d = log10(1/0)))

  ##aerobes test
  aero_idx_outcome <- aerobes_index(abx_test_df)

  testthat::expect_type(aero_idx_outcome, "double")
  testthat::expect_equal(length(aero_idx_outcome), 4)
  testthat::expect_equal(aero_idx_outcome, c(a = log10(1/0),
                                             b = log10(1/0),
                                             c = log10(1/0),
                                             d = log10(1/0)))

  ##aminoglycoside test
  aminoglycoside_idx_outcome <- aminoglycoside_index(abx_test_df)

  testthat::expect_type(aminoglycoside_idx_outcome, "double")
  testthat::expect_equal(length(aminoglycoside_idx_outcome), 4)
  testthat::expect_equal(aminoglycoside_idx_outcome, c(a = log10(1/0),
                                                       b = log10(1/0),
                                                       c = log10(1/0),
                                                       d = log10(1/0)))

})

testthat::context("Testing list of susceptible or resistant bacteria")

testthat::test_that("Testing list of susceptible or resistant bacteria", {

  ##vanco
  vanco_list <- vancomycin_list(abx_test_df)
  testthat::expect_equal(paste0(colnames(vanco_list), collapse = ", "), "SampleID, phenotype, abundance, lineage")
  testthat::expect_equal(ncol(vanco_list), 4)
  testthat::expect_equal(class(vanco_list), "data.frame")

  ##tet
  tet_list <- tetracycline_list(abx_test_df)
  testthat::expect_equal(paste0(colnames(tet_list), collapse = ", "), "SampleID, phenotype, abundance, lineage")
  testthat::expect_equal(ncol(tet_list), 4)
  testthat::expect_equal(class(tet_list), "data.frame")

  ##gram_pos
  gp_list <- gram_pos_list(abx_test_df)
  testthat::expect_equal(paste0(colnames(gp_list), collapse = ", "), "SampleID, phenotype, abundance, lineage")
  testthat::expect_equal(ncol(gp_list), 4)
  testthat::expect_equal(class(gp_list), "data.frame")

  ##gram_neg
  gn_list <- gram_neg_list(abx_test_df)
  testthat::expect_equal(paste0(colnames(gn_list), collapse = ", "), "SampleID, phenotype, abundance, lineage")
  testthat::expect_equal(ncol(gn_list), 4)
  testthat::expect_equal(class(gn_list), "data.frame")

  ##anaerobe
  anaero_list <- anaerobe_list(abx_test_df)
  testthat::expect_equal(paste0(colnames(anaero_list), collapse = ", "), "SampleID, phenotype, abundance, lineage")
  testthat::expect_equal(ncol(anaero_list), 4)
  testthat::expect_equal(class(anaero_list), "data.frame")

  ##aerobe
  aero_list <- aerobe_list(abx_test_df)
  testthat::expect_equal(paste0(colnames(aero_list), collapse = ", "), "SampleID, phenotype, abundance, lineage")
  testthat::expect_equal(ncol(aero_list), 4)
  testthat::expect_equal(class(aero_list), "data.frame")

  ##aminoglycoside
  amino_list <- aminoglycoside_list(abx_test_df)
  testthat::expect_equal(paste0(colnames(amino_list), collapse = ", "), "SampleID, phenotype, abundance, lineage")
  testthat::expect_equal(ncol(amino_list), 4)
  testthat::expect_equal(class(amino_list), "data.frame")

})


testthat::context("Testing abx_idx_df")

testthat::test_that("Testing abx_idx_df contains the required columns", {

  testthat::expect_equal(paste0(colnames(abx_idx_df), collapse = ", "), "attribute, boo, name, rank, doi")
  testthat::expect_equal(ncol(abx_idx_df), 5)
  testthat::expect_equal(class(abx_idx_df), "data.frame")

})

testthat::test_that("Testing for duplications in abx_idx_df", {

  testthat::expect_equal(sum(duplicated(abx_idx_df[, c("attribute", "name", "rank")])), 0)

})
