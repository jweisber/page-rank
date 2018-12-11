summarize_placements <- function(df, years, types) {
  df %>%
    distinct() %>%
    filter(!is.na(year) & !is.na(placement_program) & !is.na(type)) %>%
    filter(year %in% years) %>%
    filter(type %in% types) %>%
    group_by(grad_program, placement_program) %>%
    summarize(number = n())
}

square_df <- function(df, programs) {
  grid <- expand.grid(grad_program = programs, placement_program = programs, stringsAsFactors = FALSE)
  grid %>% 
    left_join(df, by = c("grad_program", "placement_program")) %>%
    mutate(number = if_else(is.na(number), 0, as.numeric(number))) %>%
    spread(key = placement_program, value = number)
}

markov_matrix <- function(df, selfies = FALSE) {
  df <- df %>% select(-1)          # remove grad_program column
  A <- data.matrix(df)             # convert to matrix
  if(!selfies) { diag(A) <- 0 }    # zero out the diagonal
  A <- t(t(A) / colSums(A))        # normalize columns
  A[is.nan(A)] <- 1 / dim(A)[1]    # set zero columns to 1/N
  return(A)
}

oneish <- function(v) {
  abs(v - 1) < 1e-15
}

page_rank <- function(A, d = .9) {
  M <- d * A + (1 - d) / dim(A)[1]

  if (!is.matrix(M)) { stop("page_rank() requires a matrix as first argument.") }
  if (dim(M)[1] != dim(M)[2]) { stop("page_rank() requires a square matrix.") }
  if (!all(oneish(colSums(M)))) { stop("page_rank() requires a column-stochastic matrix.") }
  if (sum(M <= 0)) { stop("page_rank() requires a positive matrix.") }
  
  v <- Re(eigen(M)$vectors[ , 1])
  return(v / sum(v))
}

page_rank_df <- function(df, programs, selfies = FALSE, d = .9) {
  network <- square_df(df, programs)
  A <- markov_matrix(network, selfies = selfies)
  pr <- page_rank(A, d = d)
  results <- data.frame(program = network$grad_program, page_rank = pr, stringsAsFactors = FALSE) %>% 
    arrange(desc(page_rank))
  return(results)
}
