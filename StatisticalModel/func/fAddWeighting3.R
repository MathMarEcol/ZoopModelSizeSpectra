fAddWeighting3 <- function (df, Grp) {
  # Takes a df and calculates weighting for each Type (Net vs CPR) for different Grps
  # Note that a weight of 2, for example, is equivalent to having made exactly the same observation twice. 
  # If you want to re-weight the contributions of each datum without changing the overall magnitude of the log likelihood, then you should normalize the weights (e.g. weights <- weights/mean(weights)). 
  # So if you sum the new weights, they should come to the # of observations
  Summary <- df %>% filter(Group == Grp) %>% group_by(Type) %>% 
    summarise(Num = n()) 
  df <- df %>% filter(Group == Grp) %>% 
    mutate(WtVec3 = ifelse(Type == "CPR", (Summary$Num[2] / Summary$Num[1] * X), 1 * (1 - X)), 
           WtVec3 = WtVec3 / mean(WtVec3)) %>% 
    droplevels()
}