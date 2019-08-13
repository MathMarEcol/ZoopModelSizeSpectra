fAddWeighting2 <- function (df, X) {
  # Takes a df and calculates weighting for each Tow (X% for CPR and 1-X% for Nets)
  # Takes a df and calculates weighting for each Type (Net vs CPR) for different Grps
  # Note that a weight of 2, for example, is equivalent to having made exactly the same observation twice. 
  # If you want to re-weight the contributions of each datum without changing the overall magnitude of the log likelihood, then you should normalize the weights (e.g. weights <- weights/mean(weights)). 
  # So if you sum the new weights, they should come to the # of observations
  df <- df %>%
    mutate(WtVec2 = ifelse(Type == "CPR", WtVec * X, WtVec * (100 - X)), 
           WtVec2 = WtVec2 / mean(WtVec2))
}
  