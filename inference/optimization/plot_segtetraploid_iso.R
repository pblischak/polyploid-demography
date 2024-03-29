library(tidyverse)
library(patchwork)

results_files <- list.files(path = "segtetraploid_iso/", pattern = "csv")

get_best_loglik <- function(input_file){
	df <- read_csv(paste0("segtetraploid_iso/",input_file))
	if(grepl('gbs', input_file, fixed = TRUE)){
		df$Type <- "GBS"
	} else {
		df$Type <- "WGS"
	}
	return(df[order(df$loglik, decreasing = TRUE),][1,])
}

res <- plyr::ldply(results_files, get_best_loglik)

T2_df <- data.frame(
	Type = c("GBS","GBS","WGS","WGS"),
	T2_true = c(0.25,0.5,0.25,0.5),
	value = c(0.25,0.5,0.25,0.5)
)

A <- res %>% ggplot(aes(x=Type, y=T2_est)) +
	geom_boxplot(alpha = 0.5) +
	geom_jitter(alpha = 0.8, width = 0.1) + 
	facet_wrap(.~T2_true) +
	geom_errorbar(
		data=T2_df, aes(y=NULL, ymin=value, ymax=value), size = 1.25,
		position=position_dodge(), color="blue", alpha = 0.4
	) +
	theme_bw() +
	ggtitle(
		"Segmental Allotetraploid Formation Time"
	)

dij_df <- data.frame(
	Type = c("GBS","GBS","WGS","WGS"),
	dij_true = c(0.1,0.001,0.1,0.001),
	value = c(0.1,0.001,0.1,0.001)
)

B <- res %>% ggplot(aes(x=Type, y=dij_est)) +
	geom_boxplot(alpha = 0.5) +
	geom_jitter(alpha = 0.8, width = 0.1) + 
	facet_wrap(.~dij_true) +
	geom_errorbar(
		data=dij_df, aes(y=NULL, ymin=value, ymax=value), size = 1.25,
		position=position_dodge(), color="blue", alpha = 0.4
	) +
	theme_bw() +
	ggtitle(
		"Homoeologous Exchange Rate"
	)

A + B
ggsave("segtetraploid_iso_optimization.pdf", width = 10, height = 8)
