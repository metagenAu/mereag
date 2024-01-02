
posterior_effect<-
  function(x,
           eff_col="Treatmentmetagen" ){


    x$stan_summary %>%
      as.data.frame %>%
      tibble::rownames_to_column(var='ModelTerm') %>%
      filter(grepl(eff_col,ModelTerm))->
      df1


    posteriors<- as.matrix(x)

    effect= posteriors[,colnames(posteriors) %in% eff_col] %>% as.vector()
    df1$intercept= x$stan_summary[1,1]
    df1$p_coef = sum(effect>0)/length(effect)
    df1

  }



fit_stan_NB<-
  function(X,y,the_formula,key_term,seed=123){

    stan_glm( the_formula,

           data = X,
           family = neg_binomial_2(),
           prior_intercept = normal(0, 1, autoscale = TRUE),
           prior = normal(0, 1, autoscale = TRUE),
           seed=seed) %>%
      posterior_effect(eff_col= key_term) %>%
      return()

  }




stan_by_sv<-
  function(svs,X,model_params,model_type='BEZI'){


    base = model_params$formula
    base = Reduce(paste, deparse(base))

    X[[model_params$keyterms]] %>%
      table() %>%
      data.frame() ->
      counts

    counts$Freq %>% mean() -> av_n
    counts$Freq %>% sd() -> sd_n
    key_term= model_params$keyterms
    purrr::map_df(
      svs,
      function(sv){
        model_call =
          sprintf(" fit_stan_%s(
             X=X,
             y='%s',
             the_formula = formula(paste0(sv,base)),
             key_term= key_term,
             )",model_type,sv)

        eval(parse(text=model_call))

      },.id='SV')->
      res
    res$av_n<- av_n
    res$sd_n <- sd_n

    res %>%
      return()

  }
