#' Calculate antibiotics index for Vancomycin
#'
#' @param abundance A relative abundance dataframe of bacterial taxons with samples as columns and the taxonomy lineage as the rownames, separated by `delim`
#' @param delim Delimiter for the taxonomy lineage rownames. Rownames MUST contain `Kingdom` down to `Species` level separation by the `delim` even if no information is available for missing ranks (e.g. k__Bacteria; p__Bacteroidetes; c__Bacteroidia; ; ; ; )
#'
#' @return The calculated antibiotics index each sample in the dataframe
#' @export
#'
#' @examples
#' vancomycin_index(abx_test_df)
#'
vancomycin_index <- function(abundance, delim = "; ") {
  idx <- c("gram_positive", "vancomycin")
  lineage <- row.names(abundance)
  suscept_vector <- is_susceptible(lineage, idx, delim)
  calc_index(abundance, suscept_vector)
}

#' Return taxons that are susceptible or resistant (unknown antibiotics phenotype are counted as resistant) to Vancomycin
#'
#' @param abundance A relative abundance dataframe of bacterial taxons with samples as columns and the taxonomy lineage as the rownames, separated by `delim`
#' @param delim Delimiter for the taxonomy lineage rownames. Rownames MUST contain `Kingdom` down to `Species` level separation by the `delim` even if no information is available for missing ranks (e.g. k__Bacteria; p__Bacteroidetes; c__Bacteroidia; ; ; ; )
#'
#' @return A dataframe of taxa abundances and whether they are susceptible or resistant to vancomycin
#' @export
#'
#' @examples
#' vancomycin_list(abx_test_df)
#'
vancomycin_list <- function(abundance, delim = "; ") {
  idx <- c("gram_positive", "vancomycin")
  lineage <- row.names(abundance)
  suscept_vector <- is_susceptible(lineage, idx, delim)
  get_phenotype(suscept_vector, abundance, lineage)
}

#' Calculate antibiotics index for Tetracycline
#'
#' @param abundance A relative abundance dataframe of bacterial taxons with samples as columns and the taxonomy lineage as the rownames, separated by `delim`
#' @param delim Delimiter for the taxonomy lineage rownames. Rownames MUST contain `Kingdom` down to `Species` level separation by the `delim` even if no information is available for missing ranks (e.g. k__Bacteria; p__Bacteroidetes; c__Bacteroidia; ; ; ; )
#'
#' @return The calculated antibiotics index each sample in the dataframe
#' @export
#'
#' @examples
#' tetracycline_index(abx_test_df)
#'
tetracycline_index <- function(abundance, delim = "; ") {
  idx <- "tetracycline"
  lineage <- row.names(abundance)
  suscept_vector <- is_susceptible(lineage, idx, delim)
  calc_index(abundance, suscept_vector)
}

#' Return taxons that are susceptible or resistant (unknown antibiotics phenotype are counted as resistant) to Tetracycline
#'
#' @param abundance A relative abundance dataframe of bacterial taxons with samples as columns and the taxonomy lineage as the rownames, separated by `delim`
#' @param delim Delimiter for the taxonomy lineage rownames. Rownames MUST contain `Kingdom` down to `Species` level separation by the `delim` even if no information is available for missing ranks (e.g. k__Bacteria; p__Bacteroidetes; c__Bacteroidia; ; ; ; )
#'
#' @return A dataframe of taxa abundances and whether they are susceptible or resistant to tetracycline
#' @export
#'
#' @examples
#' tetracycline_list(abx_test_df)
#'
tetracycline_list <- function(abundance, delim = "; ") {
  idx <- "tetracycline"
  lineage <- row.names(abundance)
  suscept_vector <- is_susceptible(lineage, idx, delim)
  get_phenotype(suscept_vector, abundance, lineage)
}

#' Calculate antibiotics index for Penicillin-like antibiotics
#'
#' @param abundance A relative abundance dataframe of bacterial taxons with samples as columns and the taxonomy lineage as the rownames, separated by `delim`
#' @param delim Delimiter for the taxonomy lineage rownames. Rownames MUST contain `Kingdom` down to `Species` level separation by the `delim` even if no information is available for missing ranks (e.g. k__Bacteria; p__Bacteroidetes; c__Bacteroidia; ; ; ; )
#'
#' @return The calculated antibiotics index each sample in the dataframe
#' @export
#'
#' @examples
#' penicillin_index(abx_test_df)
#'
penicillin_index <- function(abundance, delim = "; ") {
  idx <- "penicillin"
  lineage <- row.names(abundance)
  suscept_vector <- is_susceptible(lineage, idx, delim)
  calc_index(abundance, suscept_vector)
}

#' Return taxons that are susceptible or resistant (unknown antibiotics phenotype are counted as resistant) to Penicillin-like antibiotics
#'
#' @param abundance A relative abundance dataframe of bacterial taxons with samples as columns and the taxonomy lineage as the rownames, separated by `delim`
#' @param delim Delimiter for the taxonomy lineage rownames. Rownames MUST contain `Kingdom` down to `Species` level separation by the `delim` even if no information is available for missing ranks (e.g. k__Bacteria; p__Bacteroidetes; c__Bacteroidia; ; ; ; )
#'
#' @return A dataframe of taxa abundances and whether they are susceptible or resistant to penicillin-like antibiotics
#' @export
#'
#' @examples
#' penicillin_list(abx_test_df)
#'
penicillin_list <- function(abundance, delim = "; ") {
  idx <- "penicillin"
  lineage <- row.names(abundance)
  suscept_vector <- is_susceptible(lineage, idx, delim)
  get_phenotype(suscept_vector, abundance, lineage)
}

#' Calculate antibiotics index targeting gram positive bacteria such as Glycopeptides, Macrolides, Oxazolidinones, Lincosamides, and Lipopeptides aside from Vancomycin (see \code{vancomycin_index})
#'
#' @param abundance A relative abundance dataframe of bacterial taxons with samples as columns and the taxonomy lineage as the rownames, separated by `delim`
#' @param delim Delimiter for the taxonomy lineage rownames. Rownames MUST contain `Kingdom` down to `Species` level separation by the `delim` even if no information is available for missing ranks (e.g. k__Bacteria; p__Bacteroidetes; c__Bacteroidia; ; ; ; )
#'
#' @return The calculated antibiotics index each sample in the dataframe
#' @export
#'
#' @examples
#' gram_pos_index(abx_test_df)
#'
gram_pos_index <- function(abundance, delim = "; ") {
  idx <- "gram_positive"
  lineage <- row.names(abundance)
  suscept_vector <- is_susceptible(lineage, idx, delim)
  calc_index(abundance, suscept_vector)
}

#' Return taxons that are susceptible or resistant (unknown antibiotics phenotype are counted as resistant) to antibiotics targeting gram positive bacteria
#'
#' @param abundance A relative abundance dataframe of bacterial taxons with samples as columns and the taxonomy lineage as the rownames, separated by `delim`
#' @param delim Delimiter for the taxonomy lineage rownames. Rownames MUST contain `Kingdom` down to `Species` level separation by the `delim` even if no information is available for missing ranks (e.g. k__Bacteria; p__Bacteroidetes; c__Bacteroidia; ; ; ; )
#'
#' @return A dataframe of taxa abundances and whether they are susceptible or resistant to anti-gram positive antibiotics
#' @export
#'
#' @examples
#' gram_pos_list(abx_test_df)
#'
gram_pos_list <- function(abundance, delim = "; ") {
  idx <- "gram_positive"
  lineage <- row.names(abundance)
  suscept_vector <- is_susceptible(lineage, idx, delim)
  get_phenotype(suscept_vector, abundance, lineage)
}

#' Calculate antibiotics index targeting gram negatives such as Polymyxin and Aztreonam
#'
#' @param abundance A relative abundance dataframe of bacterial taxons with samples as columns and the taxonomy lineage as the rownames, separated by `delim`
#' @param delim Delimiter for the taxonomy lineage rownames. Rownames MUST contain `Kingdom` down to `Species` level separation by the `delim` even if no information is available for missing ranks (e.g. k__Bacteria; p__Bacteroidetes; c__Bacteroidia; ; ; ; )
#'
#' @return The calculated antibiotics index each sample in the dataframe
#' @export
#'
#' @examples
#' gram_neg_index(abx_test_df)
#'
gram_neg_index <- function(abundance, delim = "; ") {
  idx <- "gram_positive"
  lineage <- row.names(abundance)
  suscept_vector <- is_susceptible(lineage, idx, delim)
  calc_index(abundance, !suscept_vector)
}

#' Return taxons that are susceptible or resistant (unknown antibiotics phenotype are counted as resistant) to antibiotics targeting gram negative bacteria
#'
#' @param abundance A relative abundance dataframe of bacterial taxons with samples as columns and the taxonomy lineage as the rownames, separated by `delim`
#' @param delim Delimiter for the taxonomy lineage rownames. Rownames MUST contain `Kingdom` down to `Species` level separation by the `delim` even if no information is available for missing ranks (e.g. k__Bacteria; p__Bacteroidetes; c__Bacteroidia; ; ; ; )
#'
#' @return A dataframe of taxa abundances and whether they are susceptible or resistant to anti-gram negative antibiotics
#' @export
#'
#' @examples
#' gram_neg_list(abx_test_df)
#'
gram_neg_list <- function(abundance, delim = "; ") {
  idx <- "gram_positive"
  lineage <- row.names(abundance)
  suscept_vector <- !is_susceptible(lineage, idx, delim)
  get_phenotype(suscept_vector, abundance, lineage)
}

#' Calculate antibiotics index targeting anaerobes such as Nitroimidazole
#'
#' @param abundance A relative abundance dataframe of bacterial taxons with samples as columns and the taxonomy lineage as the rownames, separated by `delim`
#' @param delim Delimiter for the taxonomy lineage rownames. Rownames MUST contain `Kingdom` down to `Species` level separation by the `delim` even if no information is available for missing ranks (e.g. k__Bacteria; p__Bacteroidetes; c__Bacteroidia; ; ; ; )
#'
#' @return The calculated antibiotics index each sample in the dataframe
#' @export
#'
#' @examples
#' anaerobes_index(abx_test_df)
#'
anaerobes_index <- function(abundance, delim = "; ") {
  idx <- "anaerobe"
  lineage <- row.names(abundance)
  suscept_vector <- is_susceptible(lineage, idx, delim)
  calc_index(abundance, suscept_vector)
}

#' Return taxons that are susceptible or resistant (unknown antibiotics phenotype are counted as resistant) to antibiotics targeting anaerobes
#'
#' @param abundance A relative abundance dataframe of bacterial taxons with samples as columns and the taxonomy lineage as the rownames, separated by `delim`
#' @param delim Delimiter for the taxonomy lineage rownames. Rownames MUST contain `Kingdom` down to `Species` level separation by the `delim` even if no information is available for missing ranks (e.g. k__Bacteria; p__Bacteroidetes; c__Bacteroidia; ; ; ; )
#'
#' @return A dataframe of taxa abundances and whether they are susceptible or resistant to anti-anaerobic antibiotics
#' @export
#'
#' @examples
#' anaerobe_list(abx_test_df)
#'
anaerobe_list <- function(abundance, delim = "; ") {
  idx <- "anaerobe"
  lineage <- row.names(abundance)
  suscept_vector <- is_susceptible(lineage, idx, delim)
  get_phenotype(suscept_vector, abundance, lineage)
}

#' Calculate antibiotics index targeting aerobes such as Fluoroquinolone
#'
#' @param abundance A relative abundance dataframe of bacterial taxons with samples as columns and the taxonomy lineage as the rownames, separated by `delim`
#' @param delim Delimiter for the taxonomy lineage rownames. Rownames MUST contain `Kingdom` down to `Species` level separation by the `delim` even if no information is available for missing ranks (e.g. k__Bacteria; p__Bacteroidetes; c__Bacteroidia; ; ; ; )
#'
#' @return The calculated antibiotics index each sample in the dataframe
#' @export
#'
#' @examples
#' aerobes_index(abx_test_df)
#'
aerobes_index <- function(abundance, delim = "; ") {
  idx <- "aerobe"
  lineage <- row.names(abundance)
  suscept_vector <- is_susceptible(lineage, idx, delim)
  calc_index(abundance, suscept_vector)
}

#' Return taxons that are susceptible or resistant (unknown antibiotics phenotype are counted as resistant) to antibiotics targeting aerobes
#'
#' @param abundance A relative abundance dataframe of bacterial taxons with samples as columns and the taxonomy lineage as the rownames, separated by `delim`
#' @param delim Delimiter for the taxonomy lineage rownames. Rownames MUST contain `Kingdom` down to `Species` level separation by the `delim` even if no information is available for missing ranks (e.g. k__Bacteria; p__Bacteroidetes; c__Bacteroidia; ; ; ; )
#'
#' @return A dataframe of taxa abundances and whether they are susceptible or resistant to anti-aerobic antibiotics
#' @export
#'
#' @examples
#' aerobe_list(abx_test_df)
#'
aerobe_list <- function(abundance, delim = "; ") {
  idx <- "aerobe"
  lineage <- row.names(abundance)
  suscept_vector <- is_susceptible(lineage, idx, delim)
  get_phenotype(suscept_vector, abundance, lineage)
}

#' Calculate antibiotics index targeting gram-negative aerobes such as aminoglycoside
#'
#' @param abundance A relative abundance dataframe of bacterial taxons with samples as columns and the taxonomy lineage as the rownames, separated by `delim`
#' @param delim Delimiter for the taxonomy lineage rownames. Rownames MUST contain `Kingdom` down to `Species` level separation by the `delim` even if no information is available for missing ranks (e.g. k__Bacteria; p__Bacteroidetes; c__Bacteroidia; ; ; ; )
#'
#' @return The calculated antibiotics index each sample in the dataframe
#' @export
#'
#' @examples
#' aminoglycoside_index(abx_test_df)
#'
aminoglycoside_index <- function(abundance, delim = "; ") {
  gram_neg_idx <- "gram_positive"
  aerobe_idx <- "aerobe"
  aminoglycoside_idx <- "aminoglycoside"
  lineage <- row.names(abundance)
  gram_neg_vector <- !is_susceptible(lineage, gram_neg_idx, delim)
  aerobe_vector <- is_susceptible(lineage, aerobe_idx, delim)
  aminoglycoside_vector <- is_susceptible(lineage, aminoglycoside_idx, delim)
  suscept_vector <- gram_neg_vector&aerobe_vector|aminoglycoside_vector
  calc_index(abundance, suscept_vector)
}

#' Return taxons that are susceptible or resistant (unknown antibiotics phenotype are counted as resistant) to antibiotics targeting gram-negative aerobes
#'
#' @param abundance A relative abundance dataframe of bacterial taxons with samples as columns and the taxonomy lineage as the rownames, separated by `delim`
#' @param delim Delimiter for the taxonomy lineage rownames. Rownames MUST contain `Kingdom` down to `Species` level separation by the `delim` even if no information is available for missing ranks (e.g. k__Bacteria; p__Bacteroidetes; c__Bacteroidia; ; ; ; )
#'
#' @return A dataframe of taxa abundances and whether they are susceptible or resistant to antibiotics targeting gram-negative aerobes
#' @export
#'
#' @examples
#' aminoglycoside_list(abx_test_df)
#'
aminoglycoside_list <- function(abundance, delim = "; ") {
  gram_neg_idx <- "gram_positive"
  aerobe_idx <- "aerobe"
  aminoglycoside_idx <- "aminoglycoside"
  lineage <- row.names(abundance)
  gram_neg_vector <- !is_susceptible(lineage, gram_neg_idx, delim)
  aerobe_vector <- is_susceptible(lineage, aerobe_idx, delim)
  aminoglycoside_vector <- is_susceptible(lineage, aminoglycoside_idx, delim)
  suscept_vector <- gram_neg_vector&aerobe_vector|aminoglycoside_vector
  get_phenotype(suscept_vector, abundance, lineage)
}

#' Calculate the susceptiblity vector from taxonomic lineage based on the given bacterial phenotype
#'
#' @param taxa The lineage name of taxonomic ranks
#' @param idx The bacterial phenotype to calculate the susceptibility
#' @param delim The delimiter separating each taxonomic rank
#'
#' @return A vector of 0 and 1, where 0 signifies a resistant taxon or lack of information and 1 signifies a susceptible taxon
#'
is_susceptible <- function(taxa, idx, delim) {

  ##Return a character vector of taxons related to the abx of interest for pattern searching, with each element representing a taxonomic rank
  return_ranked_pattern <- function(abx_idx, TF) {
    ##Filter abx_idx_df by interested abx
    abx_df <- abx_idx_df[abx_idx_df$attribute == abx_idx, ]

    ##Grab all taxonomic classifications related to abx and save it as a vector
    pattern <- unlist(lapply(taxonomic_ranks, function(x) {
      paste0(abx_df[abx_df$rank == x & abx_df$boo == TF, "name"], collapse = "|")
    }))

    pattern[pattern == ""] <- "EMPTY"
    pattern
  }

  ##Number vector to represent the highest taxonomic rank (kingdom = 1) to the lowest taxonomic rank (species = 7)
  counting <- c(1:7)

  ##Return a vector containing the lowest taxonomic level that matches any of the patterns of the corresponding ranks in abx_idx_df
  ##Susceptible taxon vectors have positive integers while resistant taxon vectors have negative integers
  get_ranked_pheno <- function(each_row, pattern, posneg) {
    boo_list <- sapply(each_row, function (x) {
      grepl(pattern[match(x, each_row)], x)
    })

    max(boo_list*counting)*posneg
  }

  ##Initialize vector to return
  suscept_vector <- rep(0, length(taxa))

  ##If abx is tetracyciline or penicillin, assume all taxa is susceptible
  if(any(grepl("tetracycline|penicillin", idx))) {
    suscept_vector <- rep(1, length(taxa))
  }

  ##Split taxon assignments into matrix; code adapted from from qiimer package
  taxonomic_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  split_matrix <- gsub("[kpcofgs]__", "", taxa)
  split_matrix <- strsplit(as.character(split_matrix), split = delim)
  max_ranks <- max(sapply(split_matrix, length))
  split_matrix <- lapply(split_matrix, function (x) {
    fill_length <- max_ranks - length(x)
    c(x, rep("", fill_length))
  })
  split_matrix <- as.matrix(do.call(rbind, split_matrix))
  colnames(split_matrix) <- taxonomic_ranks[1:ncol(split_matrix)]

  ##Get susceptible vector for each abx in loop. Phenotypes of later abx in loop overrides earlier phenotypes
  for(each_idx in idx) {

    res_pattern <- return_ranked_pattern(each_idx, FALSE)
    sus_pattern <- return_ranked_pattern(each_idx, TRUE)

    res_numbers <- apply(split_matrix, 1, get_ranked_pheno, res_pattern, -1)
    sus_numbers <- apply(split_matrix, 1, get_ranked_pheno, sus_pattern, 1)

    if(any((abs(res_numbers) == sus_numbers) & (res_numbers != 0) & (sus_numbers != 0))) {
      dup_phenotype_idx <- which(abs(res_numbers) == sus_numbers & res_numbers != 0 & sus_numbers != 0)
      stop(simpleError(paste0("Conflicting antibiotics phenotype information regarding [",
                              paste0(taxa[dup_phenotype_idx], collapse = "] AND ["),
                              "] in abx_idx_df. Please contact tuv@chop.edu or bittingerk@chop.edu about this error")))
    }

    ##Add the susceptible and resistant taxon vectors and take the sign of the lowest taxonomic level to determine susceptibility
    new_vector <- sign(sus_numbers+res_numbers)

    ##For any taxons without information for the current looping antibiotics, use the sign of the previous antibiotics index
    zero_vector_idx <- which(new_vector == 0)
    old_suscept_vector <- suscept_vector[zero_vector_idx]
    suscept_vector <- new_vector
    suscept_vector[zero_vector_idx] <- old_suscept_vector
  }

  suscept_vector[which(suscept_vector == -1)] <- 0
  suscept_vector

}

#' Calculate the index based on the susceptible vector
#'
#' @param abundance_matrix The vector of taxonomic abundances for a sample
#' @param suscept_vector The vector of how susceptible a species is to the specified antibiotics
#'
#' @return The calculated antibiotics index for a sample
#'
calc_index <- function(abundance_matrix, suscept_vector) {

  TF_suscept_vector <- suscept_vector == 1

  apply(abundance_matrix, 2, function(x) {

    sum_suscept_taxa <- sum(x[TF_suscept_vector])

    log10((1-sum_suscept_taxa)/sum_suscept_taxa)
  })

}

#' Return a sorted dataframe of resistant and susceptible bacteria to the desired antibiotics
#'
#' @param suscept_vector The output of the `is_susceptible` function: a vector of taxons that are susceptible to the specified antibiotics
#' @param abundance A dataframe of taxonomic abundances
#' @param lineage Rowname of abundance df
#'
#' @return A sorted abundance dataframe of susceptible and resistant bacteria
#'
get_phenotype <- function(suscept_vector, abundance, lineage) {

  samples <- apply(abundance, 2, function(x) {

    sorted_abundance_suscept <- sort(x[suscept_vector == 1], index.return=TRUE, decreasing = TRUE)
    susceptibles <- data.frame(lineage = character(), abundance = numeric(), phenotype = character())
    if(length(sorted_abundance_suscept$x) != 0) {
      susceptibles <- data.frame(lineage = lineage[suscept_vector == 1][sorted_abundance_suscept$ix],
                                 abundance = sorted_abundance_suscept$x,
                                 phenotype = "susceptible")
    }
    sorted_abundance_resist <- sort(x[suscept_vector == 0], index.return=TRUE, decreasing = TRUE)
    resistances <- data.frame(lineage = character(), abundance = numeric(), phenotype = character())
    if(length(sorted_abundance_resist$x) != 0) {
      resistances <- data.frame(lineage = lineage[suscept_vector == 0][sorted_abundance_resist$ix],
                                abundance = sorted_abundance_resist$x,
                                phenotype = "resistant")
    }
    rbind(susceptibles, resistances, make.row.names = FALSE)
  })

  return_df <- do.call(rbind, samples)
  return_df$SampleID <- gsub("\\..*", "", row.names(return_df))
  return_df <- return_df[, c("SampleID", "phenotype", "abundance", "lineage")]
  row.names(return_df) <- 1:nrow(return_df)
  return_df

}
