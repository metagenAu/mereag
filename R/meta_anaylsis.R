#' Calculate Meta-Analysis Effect Size using Random Effects Model
#'
#' This function performs a meta-analysis using a random effects model. It computes the random effect size, its p-value, and the lower and upper confidence limits.
#'
#' @param data A data frame containing the meta-analysis data.
#' @param sm A character string specifying the summary measure; default is 'RR' (relative risk).
#' @return A data frame with columns: RandomEffect, p_value, LCL (Lower Confidence Limit), UCL (Upper Confidence Limit), and n_studies (number of studies included).
#' @import meta
#' @export
#' @examples
#' # Example usage
#' # calc_meta_effect(data = my_data, sm = 'OR')
calc_meta_effect <-

  function(data, sm = 'RR') {


    mgen <-
      tryCatch(
        {
          metagen(TE = Estimate,
                  seTE = `Std. Error`,
                  n.e = av_n * 2,
                  data = data,
                  sm = sm,
                  fixed = FALSE,
                  random = TRUE,
                  method.tau = "REML",
                  hakn = TRUE)
        },
        error =
          function(e)
            NA)

    if (inherits(mgen, "metagen")) {
      return(
        data.frame(
          RandomEffect = mgen$TE.random,
          p_value = mgen$pval.random,
          LCL = mgen$lower.random,
          UCL = mgen$upper.random,
          n_studies = nrow(data))
      )
    }

    data.frame(
      RandomEffect = NA,
      p_value = NA,
      LCL = NA,
      UCL = NA,
      n_studies = NA)
  }
#' Run Meta-Regression for Different Taxonomic Levels
#'
#' This function aggregates results from different studies and runs meta-regression for each unique taxonomic identifier.
#'
#' @param results A list of data frames, each representing results from a different study.
#' @param sm A character string specifying the summary measure; default is 'OR' (odds ratio).
#' @return A tibble with meta-regression results for each taxonomic identifier.
#' @export
#' @examples
#' # Example usage
#' # run_metareg(results = list_of_study_results, sm = 'RR')
run_metareg<-
  function(results,sm ='OR'){

    do.call('rbind',results) -> results_df
    column_names <- c("Phylum", "Class", "Order", "Family", "Genus")

    # Check which columns are present in the data frame
    present_columns <- column_names[column_names %in% names(results_df)]

    # Join the present columns
    if (length(present_columns) > 0) {
      results_df$tax_id <- apply(results_df[, present_columns], 1, function(x) paste(x, collapse = " "))
    } else {
      results_df$tax_id <- NA
    }
    n_genes<- unique(results_df$tax_id)
    meta_res<- vector('list',length=length(n_genes))
    names(meta_res)<- n_genes
    for( gene in n_genes){

      results_df %>%
        filter(ModelTerm=='TreatmentMetagen') %>%
        filter(tax_id %in% gene) %>%
        calc_meta_effect2(sm=sm ) -> meta1

      meta_res[[gene]]<- meta1

    }
    do.call('rbind',meta_res) %>%
      tibble::rownames_to_column(var='taxa') %>%
      tibble

  }

#' Calculate Bayesian Meta-Analysis Effect Size
#'
#' This function performs a Bayesian meta-analysis and computes the effect size, its lower and upper credible limits, and the number of studies.
#'
#' @param x A data frame with columns: mean, se (standard error), and study (identifier for each study).
#' @return A data frame with columns: RandomEffect, LCL (Lower Credible Limit), UCL (Upper Credible Limit), n_studies (number of studies included).
#' @export
#' @examples
#' # Example usage
#' # calc_bayesmeta_effect(x = my_data)
calc_bayesmeta_effect <- function(x) {


  mgen <- tryCatch({
    bayesmeta(y= x$mean,
              sigma = x$se,
              label = x$study,
              mu.prior.mean=0, mu.prior.sd=3) %>%
      .$summary %>%
      as.data.frame()
  }, error = function(e) NA)


  if (inherits(mgen, 'data.frame')) {
    return(
      data.frame(
        RandomEffect = mgen[2,2],
        LCL = mgen[5,2],
        UCL = mgen[6,2],
        n_studies = nrow(x))
    )
  }else{

    data.frame(EffectSize = NA, LCL = NA, UCL = NA, n_studies = NA,fit_type=NA)
  }
}

#' Run Bayesian Meta-Analysis for Different Taxonomic Levels
#'
#' This function aggregates results from different studies and runs a Bayesian meta-analysis for each unique taxonomic identifier.
#'
#' @param results A list of data frames, each representing results from a different study.
#' @return A tibble with Bayesian meta-analysis results for each taxonomic identifier.
#' @export
#' @import bayesmeta
#' @examples
#' # Example usage
#' # run_bayesmeta(results = list_of_study_results)
run_bayesmeta<-
  function(results){

    do.call('rbind',results) -> results_df
    column_names <- c("Phylum", "Class", "Order", "Family", "Genus")

    # Check which columns are present in the data frame
    present_columns <- column_names[column_names %in% names(results_df)]

    # Join the present columns
    if (length(present_columns) > 0) {
      results_df$tax_id <- apply(results_df[, present_columns], 1, function(x) paste(x, collapse = " "))
    } else {
      results_df$tax_id <- NA
    }
    n_genes<- unique(results_df$tax_id)
    meta_res<- vector('list',length=length(n_genes))
    names(meta_res)<- n_genes
    for( gene in n_genes){

      results_df %>%
        filter(ModelTerm=='TreatmentMetagen') %>%
        filter(tax_id %in% gene) %>%
        calc_bayesmeta_effect() -> meta1

      meta_res[[gene]]<- meta1

    }
    map_df(meta_res,~return(.x),.id='taxa')

  }

