library(tidyverse)

results_files <- list.files(path = "allotet_boots/", pattern = "csv")

get_best_loglik <- function(input_file){
	df <- read_csv(paste0("allotet_boots/",input_file))
	return(df[order(df$loglik, decreasing = TRUE),][1,])
}

res <- plyr::ldply(results_files, get_best_loglik)

q05 <- apply(res, 2, function(x) quantile(x,0.05))
print(q05)

q95 <- apply(res, 2, function(x) quantile(x,0.95))
print(q95)

avg <- apply(res, 2, mean)
stdev <- apply(res, 2, sd)

print(avg - stdev)
print(avg + stdev)
