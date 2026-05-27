test_that("DArTSet initializes correctly with SNP dataset", {
  gl <- readRDS(testthat::test_path("data", "gl_snp.RDS"))
  dartset_snp <- DArTSet$new(dataset = gl, ref=testthat::test_path("data","Pv21.fa"))
  dartset_snp$map_tags2ref(bbmap_dir = "~/Software/BBTools-39.84/")
  mta <- MTA$new(
    alleleID = "11111111-23-T/A",
    alleleSequence = "TGCAGCTTTCAACTCACCAAACATGCGACCAAACTCTTATGTCCAGAAGGGAGCACGTGTTTGGAAGCA",
    favorableAllele = "A"
  )

  expect_equal(dartset_snp$locNames, adegenet::locNames(gl))
  expect_equal(dartset_snp$indNames, adegenet::indNames(gl))
  expect_equal(dartset_snp$find_homolog(mta, cutoff = 5), c("109760619-23-T/A" = 0))
  expect_equal(dartset_snp$get_favourable_freq("109760619-23-T/A",
                                               favorable_allele = "A"), c(A=1))
  expect_equal(dartset_snp$get_favourable_freq("109760619-23-T/A",
                                               favorable_allele = "T"), c(T=0))
})

test_that("DArTSet initializes correctly with silico dataset", {
  gl <- readRDS(testthat::test_path("data", "gl_silico.RDS"))
  dartset_sil <- DArTSet$new(dataset = gl,ref=testthat::test_path("data","Pv21.fa"), dataset_type = "SIL")
  dartset_sil$map_tags2ref(bbmap_dir = "~/Software/BBTools-39.84/")
  mta <- MTA$new(
    alleleID = "11111111",
    alleleSequence = "TGCAGGATGGCAGAATGACACCCACCTCTGAACAGCTTGTACTCCAATTAGTTAATACTATTTCTTGTG",
    favorableAllele = 1
  )

  expect_equal(dartset_sil$locNames, adegenet::locNames(gl))
  expect_equal(dartset_sil$indNames, adegenet::indNames(gl))
  expect_equal(dartset_sil$find_homolog(mta, cutoff = 5), c("100454489" = 0))
  expect_equal(dartset_sil$get_favourable_freq("100454489",
                                               favorable_allele = "1"), c("1"=round(0.002747253, 6)))
  expect_equal(dartset_sil$get_favourable_freq("100454489",
                                               favorable_allele = "0"), c("0"=round(0.997252747, 6)))
})
