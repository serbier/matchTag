#' Class MTA
#'
#' Represent all metadata of a marker trait association
#'
#' @examples
#'  mta <- MTA$new(
#'    alleleID = "120515446-16-G/T",
#'    alleleSequence = "TGCAGAGGAAGAATCATGGGTTAGCTTGTTCAAATCTTCGGGTAACTAAACAAATAAAACGACTTGTGC",
#'    favorableAllele = "T"
#'  )
#'
#' @export
MTA <- R6::R6Class(
  "MTA",
  public = list(
    #' @field alleleID a string in the format of "chr-pos-allele1/allele2" for SNP MTA, and "chr-pos" for SIL MTA
    alleleID = NULL,
    #' @field alleleSequence a string of the allele sequence, with the same length as the reference genome sequence. The allele at the position of the locus tag must be either allele1 or allele2 for SNP MTA, and must be 0 or 1 for SIL MTA.
    alleleSequence = NULL,
    #' @field favorableAllele a string of the favorable allele, which must be either allele1 or allele2 for SNP MTA, and must be 0 or 1 for SIL MTA.
    favorableAllele = NULL,
    #' @field mtaType a string indicating the type of MTA, either "SNP" or "SIL"
    mtaType = NULL,
    #' @field alleles a vector of alleles for the MTA
    alleles = NULL,
    #' @field loc_tag_pos a numeric indicating the position of the locus tag in the allele sequence for SNP MTA, and NULL for SIL MTA
    loc_tag_pos = NULL,

    #' @description
    #' Initialize a new MTA object
    #'
    #' @param alleleID a string in the format of "chr-pos-allele1/allele2" for SNP MTA, and "chr-pos" for SIL MTA
    #' @param alleleSequence a string of the allele sequence, with the same length as the reference genome sequence. The allele at the position of the locus tag must be either allele1 or allele2 for SNP MTA, and must be 0 or 1 for SIL MTA.
    #' @param favorableAllele a string of the favorable allele, which must be either allele 1 or allele2 for SNP MTA, and must be 0 or 1 for SIL MTA.
    #' @return a new MTA object
    #' @export
    initialize = function(alleleID, alleleSequence, favorableAllele) {
      self$alleleID <- alleleID
      self$alleleSequence <- alleleSequence
      self$favorableAllele <- favorableAllele

      self$mtaType <- private$get_mta_type()
      self$alleles <- private$get_alleles()
      self$loc_tag_pos <- private$get_tag_loc_pos()

      private$validate_allele()
    }
  ),

  private = list(
    get_mta_type = function() {
      if (stringr::str_detect(self$alleleID, "[AGCT]/[AGCT]")) {
        return("SNP")
      } else {
        return("SIL")
      }
    },

    get_alleles = function() {
      if (self$mtaType == "SNP") {
        alleles <- self$alleleID |>
          stringr::str_extract("[AGCT]/[AGCT]$") |>
          stringr::str_split("/", simplify = T) |>
          as.vector()
      } else {
        alleles <- c("0","1")
      }
      return(alleles)
    },

    get_tag_loc_pos = function() {
      if (self$mtaType == "SNP") {
        loc_tag_pos <- self$alleleID |>
          stringr::str_split("-", simplify = T) |>
          as.vector()
        return(as.numeric(loc_tag_pos[2]))
      } else {
        return(NULL)
      }
    },

    validate_allele = function() {
      if (self$mtaType == "SNP"){
        tag_allele <- substr(self$alleleSequence,
                             self$loc_tag_pos+1,
                             self$loc_tag_pos+1)
        if (!tag_allele %in% self$alleles) {
          stop("Inconsistent allele sequence and alleles",
               paste(tag_allele, self$alleles, collapse = ", "))
        }
      } else {
        if (!self$favorableAllele %in% self$alleles) {
          stop("Favorable allele must be 0 or 1 for SIL MTA")
        }
      }
    }
  )
)
