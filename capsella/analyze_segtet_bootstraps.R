library(tidyverse)

results_files <- list.files(path = "segtet_boots/", pattern = "csv")

get_best_loglik <- function(input_file){
  df <- read_csv(paste0("segtet_boots/",input_file))
  return(df[order(df$loglik, decreasing = TRUE),][1,])
}

res <- plyr::ldply(results_files, get_best_loglik)

q025 <- apply(res, 2, function(x) quantile(x,0.025))
print(q025)

q975 <- apply(res, 2, function(x) quantile(x,0.975))
print(q975)

avg <- apply(res, 2, mean)
stdev <- apply(res, 2, sd)

print(avg - stdev)
print(avg + stdev)
