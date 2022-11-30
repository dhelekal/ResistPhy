library(CostOfResistance)
library(ggplot2)
library(reshape2)
library(patchwork)
library(hexbin)

ddf <- out$draws_df

lin1_names <- paste0("q_", c("u","t"),"[1]")
lin2_names <- paste0("q_", c("u","t"),"[2]")

lin1_df <- ddf[,lin1_names]
colnames(lin1_df) <- c("q_u","q_t")
lin1_df$lineage <- "1"

lin2_df <- ddf[,lin2_names]
colnames(lin2_df) <- c("q_u","q_t")
lin2_df$lineage <- "2"

plot_df <- rbind(lin1_df,lin2_df)

h1 <- ggplot(plot_df, aes)

