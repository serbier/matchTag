#' DArTSet class
#'
#' Represent a DArTSeq dataset of SNPs or SIL
#' @export
DArTSet <- R6::R6Class(
  "DArTSet",
  public = list(
    #'@field dataset_type A character indicating if the DArTSeq type
    dataset_type = "NULL",
    #'@field dataset A DArTSet object containing the data and metadata for the DArT analysis.
    dataset = NULL,
    #' @field locNames A character vector containing the names of the loci in the dataset.
    locNames = NULL,
    #' @field indNames A character vector containing the names of the individuals in the dataset.
    indNames = NULL,
    #' @field ref Path to reference fasta file where loci will be alligned
    ref = NULL,
    #' @field db A character vector containing the tag sequences in the dataset
    db = NULL,
    #' @field temp_dir A character string containing the path to the temporary directory.
    temp_dir = NULL,
    #' @field aligments A dataframe with mapping metadata for each tag
    aligments = NULL,

    #' @description
    #' Create a new DArTSet object.
    #' @param dataset A DArTSet object containing the data and metadata for the DArT analysis.
    #' @param dataset_type A character indicating if is "SNP" or "SIL" dataset (SNP, default).
    #' @param ref Path to a uncompressed fasta file where loci tags will be aligned.
    #' @return A DArTSet object.
    #' @export
    initialize = function(dataset, dataset_type = "SNP", ref = NULL) {
      self$dataset <- dataset
      self$dataset_type = dataset_type
      self$locNames <- adegenet::locNames(self$dataset)
      self$indNames <- adegenet::indNames(self$dataset)
      self$db <- self$build_db()
      self$temp_dir <- private$create_temp_dir()

      if(!is.null(ref)) {
        self$ref <- ref
      }
      # Useful for join mappings with tags
      self$dataset@other$loc.metrics$loc_id <- self$locNames
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
      db <- ape::as.DNAbin(mat)
      rownames(db) <- self$locNames
      return(db)
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
      freqs <- as.vector(unlist(round(gl.alf(sgl)[1,], 6)))
      names(freqs) <- as.character(imta$alleles)
      return(freqs[favorable_allele])

    },
    map_tags2ref = function(bbmap_dir, memory = "4g", use_index = T, threads = 4) {
      cp_dir <- file.path(bbmap_dir, "current/")

      # Create fasta file with tags
      private$get_tagfasta()
      in_fasta <- file.path(self$temp_dir, "tags.fa")
      out_sam <- file.path(self$temp_dir, "aligments.sam")
      idx_flag <- ifelse(use_index, "f", "t")
      args <- c(
        paste0("-Xmx", memory),
        "-cp",
        private$quote_arg(private$norm_win(cp_dir, mustWork = TRUE)),
        "align2.BBMap",
        paste0("in=",  private$quote_arg(private$norm_win(in_fasta, mustWork = TRUE))),
        paste0("out=", private$quote_arg(private$norm_win(out_sam, mustWork = FALSE))),
        paste0("ref=", private$quote_arg(private$norm_win(self$ref, mustWork = TRUE))),
        "sam=1.3",
        "ambiguous=toss",
        paste0("rebuild=", idx_flag),
        paste0("threads=", threads)
      )
      print(args)
      status <- system2("java", args)

      if (!identical(status, 0L)) {
        stop("BBDuk failed with exit status: ", status, call. = FALSE)
      }
      return(private$read_sam(min_mapq = 20))
    }

  ),
  private = list(
    create_temp_dir = function() {
      path <- tempdir()
      return(path)
    },
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
    },

    get_tagfasta = function() {
      tag_path <- file.path(self$temp_dir, "tags.fa")
      ape::write.FASTA(self$db, tag_path)
      if (file.exists(tag_path)) {
        cli::cli_inform("Tags fasta was created in {tag_path}")
      } else {
        cli::cli_abort("Failed to create tags fasta file at {tag_path}")
      }
    },

    norm_win = function(x, mustWork = TRUE) {
      normalizePath(x, winslash = "\\", mustWork = mustWork)
    },

    quote_arg = function(x) {
      shQuote(x, type = "cmd")
    },

    read_sam = function(min_mapq = 20) {
      sam_df <- read.delim(
        file.path(self$temp_dir, "aligments.sam"),
        comment.char = "@",
        header = FALSE,
        stringsAsFactors = FALSE
      )

      names(sam_df)[1:13] <- c(
        "qname", "flag", "rname", "pos", "mapq", "cigar",
        "rnext", "pnext", "tlen", "seq", "qual", "NM", "AM"
      )

      filt_map <- sam_df %>%
        dplyr::filter(mapq >= min_mapq)

      join_tab <- base::merge(self$dataset@other$loc.metrics, filt_map, by.x = "loc_id",
                              by.y = "qname", all.x = T)

      if (self$dataset_type == "SNP"){
        join_tab <- join_tab %>%
        dplyr::mutate(allele_test = stringr::str_sub(AlleleSequence, SnpPosition+1, SnpPosition+1),
               query_len = stringr::str_count(AlleleSequence))
      } else {
        join_tab <- join_tab %>%
        dplyr::mutate(SnpPosition = 1,
              allele_test = stringr::str_sub(AlleleSequence, SnpPosition+1, SnpPosition+1),
               query_len = stringr::str_count(AlleleSequence))
      }

      pred_positions <- purrr::pmap(list(join_tab$cigar,join_tab$flag,
                                         join_tab$query_len, join_tab$pos,
                                         join_tab$SnpPosition,
                                         join_tab$loc_id), private$query2ref)
      vec_positions <- purrr::map_vec(pred_positions, ~ifelse(is.null(.x), NA, .x))
      join_tab$snp_position  <-  vec_positions
      self$aligments <- join_tab
      self$dataset@chromosome <- as.factor(join_tab$rname)
      self$dataset@position <- as.integer(join_tab$snp_position)
    },

    query2ref = function(cigar, flags, query_len, ref_start, target_query_pos, alleleID){
      # Convert zero pos to 1 based pos
      target_query_pos  <- target_query_pos  + 1
      # Parse CIGAR string into lengths and operations
      cigar_tuples <- unlist(regmatches(cigar, gregexpr("\\d+[MIDNSH]", cigar)))
      # Initialize Positions
      q_pos  <- 0
      r_pos  <- ref_start - 1
      if(!is.na(flags)){
        # If tag is mapped reverse
        if(flags == 16){
          target_query_pos  <- query_len - target_query_pos + 1
          cigar_tuples  <- rev(cigar_tuples)
        }
      }
      # Loop over parsed CIGAR elements
      for (element in cigar_tuples){
        len  <- as.numeric(gsub("[MIDNSH]", "", element))
        op  <- gsub("\\d+", "", element)
        switch(op,
               M = {
                 if(q_pos + len >= target_query_pos){
                   return(r_pos + (target_query_pos - q_pos))
                 } else {
                   q_pos  <- q_pos + len
                   r_pos  <- r_pos + len
                 }
               },
               I = {
                 if(q_pos + len >= target_query_pos){
                   cli::cli_warn("Target Position located over an insertion, reporting left most pos: {alleleID}")
                   return(r_pos)
                 } else {
                   q_pos  <- q_pos + len
                 }
               },
               D = {
                 r_pos  <- r_pos + len
               },
               N = {
                 r_pos  <- r_pos + len
               },
               S = {
                 if(q_pos + len >= target_query_pos){
                   cli::cli_warn("Target Position located over an softclip, reporting left most pos {alleleID}")
                   return(r_pos)
                 } else {
                   q_pos  <- q_pos + len
                 }
               },
               H = {
                 cli::cli_warn("Found a hardclip not expected returning NA {alleleID}")
                 return(NA)
               },

               cli::cli_abort('Position not found in the alignment {alleleID}')
        )
      }
    }

  )
)
