#' Analyze microbiome trials using SQL and GAMLSS
#'
#' @param trials A list of phyloseqSparse objects representing trials
#' @param keyterms Key terms for model formula
#' @param adjterms Adjustment terms for model formula
#' @param formula Formula as string
#' @param model_type Model type, options are "BEZI", "ZIBE", etc.
#' @param new_si Optionally provide new sample information data frame
#'
#' @return A list of data frames with model results for each trial
#'
#' @examples
#' trials <- list(trial1_ps, trial2_ps)
#' results <- anaylse_trials_sql(trials, keyterms = "GroupName",
#'                               adjterms = "~ Age + BMI",
#'                               formula = "Abundance ~ GroupName",
#'                               model_type = "BEZI")
#'
#' @export

anaylse_trials_sql<-
  function(trials,keyterms,adjterms,formula,model_type, new_si=NULL,prevalence_threshold=5){

 map(trials,
        function(trial){

          sql_phyloseq_by_sample(taxa_group = 'bacteria',samples=trial) %>%
            tax_glom('Genus',NArm=FALSE) %>%
            extract_phyloseq(prevalence_threshold=prevalence_threshold) ->
            input

          svs<- colnames(input$otu_table)[-1]

          dplyr::bind_cols(
            list(input$otu_table,input$sample_data)
          ) -> X

          names(svs) = svs

          model_dict(
            keyterms =keyterms,
            adjterms=adjterms,
            formula=formula,
            data=X
            ) ->
            dict1

          gamlsss_by_sv(svs=svs,X=X,model_params=dict1,model_type=model_type ) %>%
            left_join(input$tax_table,by='SV') ->
            results
          rm(input)
          gc()
          results
        }
  ) -> results

}

#' Analyze microbiome trials using GAMLSS
#'
#' @param trials A list of phyloseqSparse objects representing trials
#' @param keyterms Key terms for model formula
#' @param adjterms Adjustment terms for model formula
#' @param formula Formula as string
#' @param model_type Model type, options are "BEZI", "ZIBE", etc.
#' @param cutoff Minimum sample size cutoff, trials below cutoff return NA
#' @param new_si Optionally provide new sample information data frame
#'
#' @return A list of data frames with model results for each trial
#'
#' @examples
#' trials <- list(trial1_ps, trial2_ps)
#' results <- anaylse_trials(trials, keyterms = "GroupName",
#'                           adjterms = "~ Age + BMI",
#'                           formula = "Abundance ~ GroupName",
#'                           model_type = "BEZI")
#'
#' @export

anaylse_trials<-
  function(trials,keyterms,adjterms,formula,model_type='BEZI',cutoff =7,prevalence_threshold=5,new_si=NULL){

    TSS<- TRUE
    if(model_type %in% c('NB','ZINBI')){
      TSS<- FALSE
    }


    map(trials,
        function(trial){

          if(phyloseqSparse::nsamples(trial)> cutoff){

          trial  %>%
            tax_glom('Genus',NArm=FALSE) %>%
            extract_phyloseq(TSS=TSS,new_si=new_si,prevalence_threshold=prevalence_threshold) ->
            input

          svs<- colnames((input$otu_table)[-1])

         # input$otu_table[is.na(input$otu_table)]<- 0


          dplyr::bind_cols(
            list(input$otu_table,input$sample_data)
          ) -> X

          print(dim(X))


          X = na.omit(X[,c(keyterms,adjterms,svs,'MetagenNumber')])

          names(svs) = svs

          model_dict(
            keyterms =keyterms,
            adjterms=adjterms,
            formula=formula,
            data=X
          ) ->
            dict1
          print('dictionary created')

          gamlsss_by_sv(svs=svs,X=X,model_params=dict1,model_type=model_type) %>%
            left_join(input$tax_table,by='SV') ->
            results
          print('complete')
          rm(input)
          gc()
          results
          }
          else{
            NA
          }
        }
    ) -> results

  }

#' Analyze microbiome trials using rstanarm and GAMLSS
#'
#' @param trials A list of phyloseqSparse objects representing trials
#' @param keyterms Key terms for model formula
#' @param adjterms Adjustment terms for model formula
#' @param formula Formula as string
#' @param model_type Model type, options are "NB", "ZINB", etc.
#' @param cutoff Minimum sample size cutoff, trials below cutoff return NA
#'
#' @return A list of data frames with model results for each trial
#'
#' @examples
#' trials <- list(trial1_ps, trial2_ps)
#' results <- anaylse_trials_stan(trials, keyterms = "GroupName",
#'                                adjterms = "~ Age + BMI",
#'                                formula = "Abundance ~ GroupName",
#'                                model_type = "ZINB")
#'
#' @export
anaylse_trials_stan<-
  function(trials,keyterms,adjterms,formula,model_type='BEZI',cutoff =7,prevalence_threshold=5){

    TSS<- TRUE
    if(model_type %in% c('NB','ZINBI')){
      TSS<- FALSE
    }


    map(trials,
        function(trial){

          if(phyloseqSparse::nsamples(trial)> cutoff){

            trial  %>%
              tax_glom('Genus',NArm=FALSE) %>%
              extract_phyloseq(TSS=TSS,prevalence_threshold=prevalence_threshold) ->
              input

            svs<- colnames(input$otu_table)[-1]

            dplyr::bind_cols(
              list(input$otu_table,input$sample_data)
            ) -> X

            names(svs) = svs

            model_dict(
              keyterms =keyterms,
              adjterms=adjterms,
              formula=formula,
              data=X

            ) ->
              dict1
            print('dictionary created')

            stan_by_sv(svs=svs,X=X,model_params=dict1,model_type=model_type) %>%
              left_join(input$tax_table,by='SV') ->
              results
            rm(input)
            gc()
            results
          }
          else{
            NA
          }
        }
    ) -> results

  }
