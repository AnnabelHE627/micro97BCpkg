usethis::use_package(package="TwoSampleMR",type="Imports")
usethis::use_package(package="mr.raps",type="Imports")
usethis::use_package(package="tidyverse",type="Depends")

#' Format data for MR analysis
#'
#' @param x Output from \code{\link[TwoSampleMR]{format_data}}. Must have a SNP name column (SNP), SNP chromosome column (chr_name), SNP position column (chrom_start). If id.exposure or pval.exposure not present they will be generated.
#'
#' @return A data frame with the same columns as x, but with only the SNPs that are insignificant at p<1E-5 and in LD with each other at r2<0.001 removed.
#' @export
#' @importFrom TwoSampleMR format_data clump_data
#'
format_ins_select<-function(x){
  filter(x,pval<1E-5) %>%
    format_data(type="exposure",
                snp_col="SNP",beta_col="beta",phenotype_col="exposure",
                se_col="se",effect_allele_col="effect_allele",other_allele_col="other_allele",
                eaf_col = "eaf",pval_col = "pval",samplesize_col = "N")%>%
    clump_data(clump_r2 = 0.001,pop="EAS")
}


#' MR_custom
#'
#' @param exp_dat exposure data
#' @param out_dat outcome data
#' @param outcome_gwas_name name of outcome GWAS
#'
#' @return A list of MR results
#' @export
#' @importFrom TwoSampleMR harmonise_data mr mr_pleiotropy_test mr_heterogeneity mr_method_list
#' @importFrom mr.raps mr.raps.simple
#'
MR_custom<-function(exp_dat,out_dat,outcome_gwas_name){
  dat<-result_MR<-result_pleiotropic<-result_heterogeneity<-result_mrraps<-list(NA)
  for (i in 1:length(exp_dat)) {
    dat[[i]]<-harmonise_data(exposure_dat = exp_dat[[i]],
                             outcome_dat = out_dat[[i]],
                             action = 2)
  }

  for (i in 1:length(exp_dat)) {
    result_MR[[i]]<-mr(dat[[i]])
    result_pleiotropic[[i]]<-mr_pleiotropy_test(dat[[i]])
    result_heterogeneity[[i]]<-mr_heterogeneity(dat[[i]])
    result_mrraps[[i]]<-mr.raps.simple(b_exp = dat[[i]]$beta.exposure,
                                       b_out = dat[[i]]$beta.outcome,
                                       se_exp = dat[[i]]$se.exposure,
                                       se_out = dat[[i]]$se.outcome)%>%t()
  }
  z<-list(NA)
  z$result_MR<-result_MR
  z$result_pleiotropic<-result_pleiotropic
  z$result_heterogeneity<-result_heterogeneity
  z$result_mrraps<-result_mrraps
  z
}
