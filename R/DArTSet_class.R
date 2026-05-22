#' DArTSet class
#'
#' Represent a DArTSeq dataset of SNPs or SIL
#' @export
DArTSet <- R6::R6Class(
  "DArTSet",
  public = list(
    #'@field dataset A DArTSet object containing the data and metadata for the DArT analysis.
    dataset = NULL,
    #' @field locNames A character vector containing the names of the loci in the dataset.
    locNames = NULL,
    #' @field indNames A character vector containing the names of the individuals in the dataset.
    indNames = NULL,
    #' @field db A character vector containing the tag sequences in the dataset
    db = NULL,


    #' @description
    #' Create a new DArTSet object.
    #' @param dataset A DArTSet object containing the data and metadata for the DArT analysis.
    #' @return A new DArTSet object.
    #' @export
    initialize = function(dataset) {
      self$dataset <- dataset
      self$locNames <- adegenet::locNames(self$dataset)
      self$indNames <- adegenet::indNames(self$dataset)
      self$db <- self$build_db()
    },

    #' @description
    #' Build a DNA sequence database from marker allele sequences.
    #'
    #' Converts `dataset@other$loc.metrics$AlleleSequence` into a DNA binary
    #' matrix (`DNAbin`) indexed by locus names.
    #'
    #' @return A `DNAbin` object containing one aligned sequence per locus.
    build_db = function() {
      mat <- do.call(rbind,
                     strsplit(as.vector(self$dataset@other$loc.metrics$AlleleSequence), ""))
      rownames(mat) <- self$lodNames
      return(ape::as.DNAbin(mat))
    },

    #' @description
    #' Find homologous loci in the dataset for a given MTA.
    #'
    #' Computes Hamming distance from the MTA allele sequence to all loci in the
    #' internal database. For SNP MTAs, matches are additionally filtered by SNP
    #' position.
    #'
    #' @param mta An `MTA` object containing the query allele metadata.
    #' @param cutoff Maximum Hamming distance accepted to call a homolog.
    #'
    #' @return A named numeric vector of distances for matches, or `NA` when no
    #'   matching locus is found.
    find_homolog = function(mta, cutoff = 5) {

      d <- private$hamming_dnabin(mta$alleleSequence)
      match_idx <- which(d <= cutoff)
      if (length(match_idx) > 0){
        matches <- d[match_idx]
        hom_names <- names(matches)
        if (mta$mtaType == "SNP") {
          possible_snp_pos <- as.numeric(stringr::str_split(hom_names, "-", simplify = T)[,2])
          pos_match_idx <- which(possible_snp_pos %in% mta$loc_tag_pos)
          match_names <- hom_names[pos_match_idx]
          if (length(match_names) > 0){
            return(d[match_names])
          } else { return(NA) }
        } else {
          return(d[match_idx])
        }
      } else {
        return(NA)
      }
    },

    #' @description
    #' Get the allele frequency of the favorable allele for one locus.
    #'
    #' Builds an `MTA` object for the requested locus, validates the favorable
    #' allele, and returns its frequency from genotype means.
    #'
    #' @param alleleID Marker identifier in `locNames`.
    #' @param favorable_allele Favorable allele state to retrieve (`0`, `1`, or
    #'   SNP allele code depending on marker type).
    #'
    #' @return Numeric scalar with the favorable allele frequency.
    get_favourable_freq = function(alleleID, favorable_allele) {
      idx <- private$get_loc_idx(alleleID)
      sgl <- self$dataset[,idx]

      imta <- MTA$new(alleleID = alleleID,
                      alleleSequence = sgl@other$loc.metrics$AlleleSequence,
                      favorableAllele = favorable_allele)
      if (!favorable_allele %in% imta$alleles){
        stop(paste("Favorable allele", favorable_allele,
                   "not found in the dataset for alleleID", alleleID))
      }
      freq_alt <- as.vector(adegenet::glMean(sgl))
      freq_ref <- 1 - freq_alt
      freqs <- c(freq_ref, freq_alt)
      names(freqs) <- imta$alleles
      return(freqs[favorable_allele])

    }

  ),
  private = list(
    get_loc_idx = function(alleleID) {
      idx <- which(self$locNames == alleleID)
      if (length(idx) == 0) {
        warning(paste("AlleleID", alleleID, "not found in the dataset."))
        return(NA)
      } else {
        return(idx)
      }
    },
    hamming_dnabin = function(query) {
      qnn <- ape::as.DNAbin(unlist(strsplit(query, "")))
      qnn_raw <- as.raw(qnn)
      d <- rowSums(self$db != matrix(qnn_raw,
                            nrow = nrow(self$db),
                            ncol = ncol(self$db),
                            byrow = TRUE))
      names(d) <- self$locNames
      return(sort(d, decreasing = T))
    }
  )
)
