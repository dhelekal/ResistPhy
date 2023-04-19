library(ggplot2)
library(patchwork)
p1 <- plot_epi_dynamic(out, 2) +
    scale_x_continuous(breaks = seq(t0, tmax, by =2)) +
    ggtitle("Resistant Lineage 1")
p1 <- p1 + (p1$theme %+replace% theme(legend.position="none",
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank()))

p2 <- plot_epi_dynamic(out, 3) +
    scale_x_continuous(breaks = seq(t0, tmax, by =2)) +
    ggtitle("Resistant Lineage 2")
p2 <- p2 + (p2$theme %+replace% theme(axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank()))

y1_min <- min(layer_scales(p1)$y$range$range, layer_scales(p2)$y$range$range)
y1_max <- max(layer_scales(p1)$y$range$range, layer_scales(p2)$y$range$range)

p1 <- p1+ylim(y1_min,y1_max)
p2 <- p2+ylim(y1_min,y1_max)


p3 <- plot_rt(out, 2) +
    scale_x_continuous(breaks = seq(t0, tmax, by =2))
p3 <- p3 + (p3$theme %+replace% theme(legend.position="none",plot.title=element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))) 

p4 <- plot_rt(out, 3) +
    scale_x_continuous(breaks = seq(t0, tmax, by =2))
p4 <- p4 + (p4$theme %+replace% theme(axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    plot.title=element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))) 

fig3 <- (p1 | p2)/(p3 | p4)

