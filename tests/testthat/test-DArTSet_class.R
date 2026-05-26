test_that("DArTSet initializes correctly with SNP dataset", {
  gl <- readRDS(testthat::test_path("data", "gl_snp.RDS"))
  dartset_snp <- DArTSet$new(dataset = gl, ref="./tests/Pv21.fa.gz")
  mta <- MTA$new(
    alleleID = "109760506-56-T/C",
    alleleSequence = "TGCAGAAAACGAGAAACCCAGAACCAAAAACCAAACCAGAGCCAGATCCATACCCGTTCTGGTGGTGGC",
    favorableAllele = "T"
  )

  expect_equal(dartset_snp$locNames, adegenet::locNames(gl))
  expect_equal(dartset_snp$indNames, adegenet::indNames(gl))
  expect_equal(dartset_snp$find_homolog(mta, cutoff = 5), c("109760506-56-T/C" = 0))
  expect_equal(dartset_snp$get_favourable_freq("109760506-56-T/C",
                                               favorable_allele = "C"), c(C=1))
  expect_equal(dartset_snp$get_favourable_freq("109760506-56-T/C",
                                               favorable_allele = "T"), c(T=0))
})

test_that("DArTSet initializes correctly with silico dataset", {
  gl <- readRDS(testthat::test_path("data", "gl_silico.RDS"))
  dartset_sil <- DArTSet$new(dataset = gl)
  mta <- MTA$new(
    alleleID = "111111111",
    alleleSequence = "TGCAGCCACATTTTGCTGGGTCAATGTTCAACATGGGGGGCATCGTTTATGATTTTTTACAGATCGGAA",
    favorableAllele = 1
  )

  expect_equal(dartset_sil$locNames, adegenet::locNames(gl))
  expect_equal(dartset_sil$indNames, adegenet::indNames(gl))
  expect_equal(dartset_sil$find_homolog(mta, cutoff = 5), c("121779781" = 0))
  expect_equal(dartset_sil$get_favourable_freq("121779781",
                                               favorable_allele = "1"), c("1"=0))
  expect_equal(dartset_sil$get_favourable_freq("121779781",
                                               favorable_allele = "0"), c("0"=1))
})
