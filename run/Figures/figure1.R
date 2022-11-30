library(ggplot2)
library(latex2exp)
library(viridis)
library(hexbin)
library(patchwork)

p1 <- plot_epi_dynamic(out, 1, n_draws=20) + geom_line(data=df,aes(x=t, y=I_2), color="red") + ggtitle("Susceptible Lineage")
p1 <- p1 + (p1$theme %+replace% theme(legend.position="none",
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank()))
p2 <- plot_epi_dynamic(out, 2, n_draws=20) + geom_line(data=df,aes(x=t, y=I_3), color="red") + ggtitle("Resistant Lineage")
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

p3 <- plot_rt(out, 1, n_draws=20) + geom_line(data=df,aes(x=t, y=Rt_2), color="red")
p3 <- p3 + (p3$theme %+replace% theme(legend.position="none",plot.title=element_blank()))
p4 <- plot_rt(out, 2, n_draws=20) + geom_line(data=df,aes(x=t, y=Rt_3), color="red") 
p4 <- p4 + (p4$theme %+replace% theme(axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    plot.title=element_blank()))

y2_min <- min(layer_scales(p3)$y$range$range, layer_scales(p4)$y$range$range)
y2_max <- max(layer_scales(p3)$y$range$range, layer_scales(p4)$y$range$range)

p3 <- p3+ylim(y2_min,y2_max)
p4 <- p4+ylim(y2_min,y2_max)

fig1 <- (p1 | p2) / (p3 | p4) 