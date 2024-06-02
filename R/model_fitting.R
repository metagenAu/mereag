
#' Validate Column
#'
#' Checks if a specified column exists in the dataframe and has data.
#'
#' @param df Dataframe
#' @param column_name Column name to validate
#' @import gamlss
#' @import rstanarm
#' @import dplyr
#' @import tibble
#' @import magrittr
#' @return Logical. TRUE if column exists with data, FALSE otherwise.
#' @export

validate_column <- function(df, column_name) {

  # Check if column exists
  if(!column_name %in% names(df)) {
    warning(paste0(column_name, " does not exist in the dataframe!"))
    return(FALSE)
  }

  # Check column has values
  if(all(is.na(df[[column_name]]))) {
    warning(paste0(column_name, " has all NA values!"))
  }

  # Column checks passed
  message(paste0(column_name, " exists in the dataframe and has data."))
  return(TRUE)

}


#' Create Model Dictionary
#'
#' Validates model terms and creates model dictionary.
#'
#' @param keyterms Key term
#' @param formula Formula string
#' @param data Model data frame
#' @param adjterms Adjustment terms
#' @return List containing model terms
#' @export
model_dict <- function(keyterms, formula, data, adjterms = NULL) {
  assert_that(length(keyterms) == 1, msg = 'Key Term Arg Missing')

  key_val <- validate_column(data, keyterms)
  assert_that(key_val == TRUE, msg = 'Key Term Missing')

  # Check if the formula contains '~'

  if (!is.null(adjterms)) {
    sapply(adjterms, function(x) validate_column(data, x)) -> adj_val
    assert_that(sum(adj_val) == length(adjterms), msg = 'Adj Term Missing from dataset')
    sapply(adjterms, function(x) grepl(x, formula)) -> adj_val_form
    assert_that(sum(adj_val_form) == length(adjterms), msg = 'One or more adj terms not in formula')

    # Define a function to check if a column is a character or factor with the required criteria
    check_column <- function(column) {
      if (is.character(column) || is.factor(column)) {
        unique_values <- unique(column)
        if (length(unique_values) >= 2) {
          counts <- table(column)
          if (all(counts >= 3)) {
            return(TRUE)
          }
        }
        return(FALSE)
      }
      return(TRUE)
    }

    # Function to check if terms in the formula are valid
    check_term <- function(term, valid_columns) {
      if (term %in% valid_columns) {
        return(TRUE)
      } else if (grepl("\\*", term)) {  # Check for interaction terms using *
        interaction_terms <- unlist(strsplit(term, "\\*"))
        return(all(interaction_terms %in% valid_columns))
      }
      return(FALSE)
    }

    # Get the terms from the original formula
    original_terms <- all.vars(as.formula(formula))

    # Filter the columns in adjterms based on the criteria
    valid_columns <- names(data)[sapply(data[adjterms], check_column)]

    # Get the terms and interactions from the formula
    terms_obj <- terms(as.formula(formula))
    term_labels <- attr(terms_obj, "term.labels")

    # Check each term and preserve valid ones including interactions
    valid_terms <- sapply(term_labels, check_term, valid_columns = valid_columns)

    # Create the filtered formula
    valid_term_labels <- term_labels[valid_terms]
    predictors <- paste(valid_term_labels, collapse = " + ")
    if(nchar(predictors)>0){
      predictors <- paste0(keyterms ,'+' , predictors)
      print(predictors)
    }else{
      predictors<- keyterms

    }
    updated_formula <- as.formula(paste("~", predictors,'+offset(seq)'))


    list(keyterms = keyterms,
         adjterms = adjterms,
         formula = updated_formula) %>%
      return()

  } else {
    warning('No Adj term detected. Check your data set is correct.')

    list(keyterms = keyterms,
         adjterms = NA,
         formula = formula) %>%
      return()
  }
}


#' Fit GAMLSS Model By Species
#'
#' Fits GAMLSS to each species in a phylogenetic sequence.
#'
#' @param svs Species identifiers
#' @param X Model data frame
#' @param model_params Model dictionary
#' @param model_type Model type
#' @import gamlss
#' @import gamlss.dist
#' @return Data frame with model results for each species
#' @export
gamlsss_by_sv<-
  function(svs,X,model_params,model_type='BEZI'){


    base = model_params$formula
    base = Reduce(paste, deparse(base))

     X[[model_params$keyterms]] %>%
      table() %>%
      data.frame() ->
      counts

     counts$Freq %>% mean() -> av_n
     counts$Freq %>% sd() -> sd_n

    purrr::map_df(
      svs,
         function(sv){
           model_call =
             sprintf(" fit_gamlsss_%s(
             X=X,
             y='%s',
             the_formula = formula(paste0(sv,base))
             )",model_type,sv)

           eval(parse(text=model_call))

         },.id='SV')->
      res
    res$av_n<- av_n
    res$sd_n <- sd_n

    res %>%
          return()

  }

#' Fit BEZI GAMLSS Model
#'
#' Fits a Beta Zero-Inflated GAMLSS model using the BEZI distribution.
#' Performs stepwise model selection and extracts key model summary statistics.
#'
#' @param X Data frame containing model matrix
#' @param y Response variable name
#' @param the_formula Formula for gamlss model
#' @param key_term Name of key variable of interest
#' @param robust Logical, use robust standard errors
#'
#' @return Tibble with model summary statistics if model converges, NA otherwise
#'
#' @export
fit_gamlsss_BEZI<-
  function(X,y,the_formula,key_term,robust=FALSE){

    gamlss::gamlss(
      the_formula,
      family = BEZI(),
      data =X,
      trace = FALSE ,
      control = gamlss.control(n.cyc = 200)
    )  ->
      mod1

    stepGAIC(mod1, what="sigma", scope=formula(mod1),direction='both') -> mod2
    stepGAIC(mod2, what="nu", scope=formula(mod2),direction='both') -> gamlss_results


    if(gamlss_results$converged){

      gamlss_results %>%
        summary(robust=robust) %>%
        as_tibble(rownames='ModelTerm')-> output

      yhat = predict(gamlss_results)
      cor1 = cor(X[[y]],yhat,use='complete.obs')
      comps = c('','sigma-','nu-')
      idx =  grepl('Intercept' ,output$ModelTerm) %>% which()
      output$ModelTerm[idx] <- paste0(comps,output$ModelTerm[idx])
      output$cor<- cor1

      output  %>%
        return()


    }else{
      return(NA)
    }

  }


#' Fit Negative Binomial GAMLSS Model
#'
#' @inheritParams fit_gamlsss_BEZI
#' @export
fit_gamlsss_NB<-
  function(X,y,the_formula,key_term,robust=FALSE){

    gamlss::gamlss(
      the_formula,
      family = NBI(),
      data =X,
      trace = FALSE ,
      control = gamlss.control(n.cyc = 200)
    )  ->
      mod1

    stepGAIC(mod1, what="sigma", scope=formula(mod1),direction='both') -> gamlss_results


    if(gamlss_results$converged){

      gamlss_results %>%
        summary(robust=robust) %>%
        as_tibble(rownames='ModelTerm')-> output

      yhat = predict(gamlss_results)
      cor1 = cor(X[[y]],yhat,use='complete.obs')
      comps = c('','sigma-')
      idx =  grepl('Intercept' ,output$ModelTerm) %>% which()
      output$ModelTerm[idx] <- paste0(comps,output$ModelTerm[idx])
      output$cor<- cor1

      output  %>%
        return()


    }else{
      return(NA)
    }

  }

#' Fit Zero-Inflated Negative Binomial GAMLSS Model
#'
#' @inheritParams fit_gamlsss_BEZI
#' @export
fit_gamlsss_ZINBI<-
  function(X,y,the_formula,key_term,robust=FALSE){

    gamlss::gamlss(
      the_formula,
      family = ZINBI(),
      data =X,
      trace = FALSE ,
      control = gamlss.control(n.cyc = 200)
    ) ->
      mod1

    stepGAIC(mod1, what="sigma", scope=formula(mod1),direction='both') -> gamlss_results


    if(gamlss_results$converged){

      gamlss_results %>%
        summary(robust=robust) %>%
        as_tibble(rownames='ModelTerm')-> output

      yhat = predict(gamlss_results)
      cor1 = cor(X[[y]],yhat,use='complete.obs')
      comps = c('','sigma-')
      idx =  grepl('Intercept' ,output$ModelTerm) %>% which()
      output$ModelTerm[idx] <- paste0(comps,output$ModelTerm[idx])
      output$cor<- cor1

      output  %>%
        return()


    }else{
      return(NA)
    }

  }
#' Fit Beta Inflated GAMLSS Model
#'
#' @inheritParams fit_gamlsss_BEZI
#' @export
fit_gamlsss_BEINF0<-
  function(X,y,the_formula,key_term,robust=FALSE){

    gamlss::gamlss(
      the_formula,
      family = BEINF0(),
      data =X,
      trace = FALSE ,
      control = gamlss.control(n.cyc = 200)
    ) ->
      mod1

    stepGAIC(mod1, what="sigma", scope=formula(mod1),direction='both') -> mod2
    stepGAIC(mod2, what="nu", scope=formula(mod2),direction='both') -> gamlss_results

    if(gamlss_results$converged){

      gamlss_results %>%
        summary(robust=robust) %>%
        as_tibble(rownames='ModelTerm')-> output

      yhat = predict(gamlss_results)
      cor1 = cor(X[[y]],yhat,use='complete.obs')
      comps = c('','sigma-','nu-')
      idx =  grepl('Intercept' ,output$ModelTerm) %>% which()
      output$ModelTerm[idx] <- paste0(comps,output$ModelTerm[idx])
      output$cor<- cor1

      output  %>%
        return()


    }else{
      return(NA)
    }

  }

