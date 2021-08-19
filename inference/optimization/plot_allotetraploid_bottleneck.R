library(tidyverse)
library(patchwork)

results_files <- list.files(path = "allotetraploid_bottleneck/", pattern = "csv")

get_best_loglik <- function(input_file){
	df <- read_csv(paste0("allotetraploid_bottleneck/",input_file))
	if(grepl('gbs', input_file, fixed = TRUE)){
		df$Type <- "GBS"
	} else {
		df$Type <- "WGS"
	}
	return(df[order(df$loglik, decreasing = TRUE),][1,])
}

res <- plyr::ldply(results_files, get_best_loglik)

A <- res %>% ggplot(aes(x=Type, y=nuBot_est)) +
	geom_boxplot(alpha = 0.5) +
	geom_jitter(alpha = 0.8, width = 0.1) + 
	facet_grid(nuBot_true~T2_true) +
	theme_classic() +
	ggtitle(
		"Effective Population Size"
	)

lines_df <- data.frame(
	Type = c("GBS","GBS","WGS","WGS"),
	T2_true = c(0.5,1.5,0.5,1.5),
	value = c(0.5,1.5,0.5,1.5)
)

B <- res %>% ggplot(aes(x=Type, y=T2_est)) +
	geom_boxplot(alpha = 0.5) +
	geom_jitter(alpha = 0.8, width = 0.1) + 
	facet_wrap(.~T2_true) +
	geom_errorbar(
		data=lines_df, aes(y=NULL, ymin=value, ymax=value), size = 1.25,
		position=position_dodge(), color="blue", alpha = 0.4
	) +
	theme_classic() +
	ggtitle(
		"Divergence Time"
	)

A + B
ggsave("allotetraploid_iso_optimization.pdf", width = 10, height = 8)
