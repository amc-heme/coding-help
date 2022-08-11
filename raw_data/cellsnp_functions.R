## This file contains draft functions to:
##  - read cellsnp-style VCFs (row metadata)
##  - read cellsnp output to a consistently structured list that contains a set of matrices and metadata.

library(tidyverse)
library(Matrix)
library(Seurat)

#' Read a cellsnp VCF file.
#'
#' @param vcf vcf filename
#'
#' @family read functions
#' @return `tibble`
#'
#' @export
read_vcf_cellsnp <- function(vcf) {
  res <- suppressMessages(readr::read_tsv(vcf, comment = "##"))
  colnames(res) <- stringr::str_replace(colnames(res), "^#", "")
  
  res <- separate_rows(res, INFO, sep = ";") %>%
    separate(
      INFO,
      into = c("key", "value"),
      sep = "=",
      convert = TRUE
    ) %>%
    pivot_wider(names_from = "key", values_from = "value")
  
  res
}

read_cellsnp <- function(path, prefix = NULL, return_ALT = TRUE, return_REF = TRUE, return_OTH = TRUE, depth_incl_OTH = FALSE) {
  
  add_labels <- function(mat) {
    colnames(mat) <- barcodes
    row.names(mat) <- vcf$rowname_id
    mat
  }
  
  barcodes <- read_lines(file.path(path, "cellSNP.samples.tsv"))
  barcodes <- paste0(prefix, barcodes)
  
  vcf <-
    read_vcf_cellsnp(file.path(path, "cellSNP.base.vcf.gz")) %>% unite(rowname_id, CHROM, POS, REF, ALT, sep = "_", remove = FALSE)
  
  AD_mat <-
    readMM(file.path(path, "cellSNP.tag.AD.mtx")) %>% add_labels()
  DP_mat <-
    readMM(file.path(path, "cellSNP.tag.DP.mtx")) %>% add_labels()
  OTH_mat <-
    readMM(file.path(path, "cellSNP.tag.OTH.mtx")) %>% add_labels()
  REF_mat <- DP_mat - AD_mat

  if (depth_incl_OTH == TRUE) {DP_mat <- DP_mat + OTH_mat}
  
  ref <- REF_mat
  alt <- AD_mat
  
  ref[ref > 0] <- 1
  alt[alt > 0] <- 2
  
  genotype <- ref + alt
  
  out <- list(genotype = genotype,
              depth = DP_mat,
              variant_meta = vcf)
  
  if (return_ALT == TRUE) {out$alt <- AD_mat}
  if (return_REF == TRUE) {out$ref <- REF_mat}
  if (return_OTH == TRUE) {out$other <- OTH_mat}
  
  out
  
}
