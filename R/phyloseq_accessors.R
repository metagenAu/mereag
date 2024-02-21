
#' Filter and extract data from phyloseq
#'
#' @param physeq A phyloseq object
#' @param prevalence_threshold Prevalence threshold (proportion of samples) for filtering OTUs
#'
#' @return A named list with 3 elements:
#' \itemize{
#'   \item{\code{otu_table}}{OTU table as tibble with OTU_ID column}
#'   \item{\code{tax_table}}{Taxonomy table as tibble with Taxonomy_ID column}
#'   \item{\code{sample_data}}{Sample data as tibble}
#' }
#'
#' @examples
#' filtered_data <- extract_phyloseq(physeq, prevalence_threshold = 0.05, abundance_threshold = 0.001)
#'
#' @export
#' @import phyloseqSparse
#' @import dplyr
#' @import tidyr
extract_phyloseq <-

  function(physeq, TSS=TRUE,prevalence_threshold = 5,new_si=NULL) {

  #  sample_data(physeq)$seq = log(sample_sums(physeq))

    if(!is.null(new_si)){

      sample_data(physeq)<- new_si

    }

    physeq %>%
      transform_sample_counts(
        function(x){
          if(TSS){
            x / sum(x)
            }else{
              x
              }
          } ) %>%
      filter_taxa(function(x) sum(x>0) > prevalence_threshold , TRUE) %>%
      otu_table() %>%
      as.matrix() %>%
      as_tibble(rownames="sampleID") -> otu_df

    tax_df<- NA

    try(physeq %>%
      tax_table() %>%
      as.data.frame()  %>%
      as_tibble(rownames="SV") -> tax_df)

    physeq %>%
      sample_data() %>%
      as.data.frame %>%
      as_tibble() -> sample_df

    if(TSS==FALSE){
      sample_df$seq<- log(rowSums(otu_df[,-1]))
    }

    return(list(
      otu_table = otu_df,
      tax_table = tax_df,
      sample_data = sample_df
    ))

  }
#' Create Identifier Column
#'
#' Creates a character column that concatenates selected columns with underscores.
#' Useful for creating identifiers.
#'
#' @param tx Data frame
#' @param level Column name to concatenate up to
#'
#' @return The original data frame with the new identifier column added
#'
#' @examples
#' mydf <- data.frame(Col1 = letters[1:3],
#'                    Col2 = 1:3,
#'                    Col3 = c("Group1", "Group2", "Group1"))
#' newdf <- create_id(mydf, "Col3")
#'
#' @export
create_id<- function(tx,level){

  cols = colnames(tx)

  idx= which(cols==level)

  cols_to_keep = cols[1:idx]

  unite(tx, cols_to_keep, sep = "_", remove = TRUE, na.rm = FALSE)


}
